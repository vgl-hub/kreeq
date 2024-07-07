#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <functional>
#include <list>
#include <iomanip>
#include <stdio.h>
#include <chrono>
#include <array>
#include <atomic>
#include <future>
#include <cstdio>
#include <deque>
#include <limits>

#include "parallel-hashmap/phmap.h"
#include "parallel-hashmap/phmap_dump.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"
#include "output.h"
#include "string-graph.h"

#include "input.h"
#include "kmer.h"
#include "kreeq.h"
#include "fibonacci-heap.h"

#include "variants.h"

void DBG::correctSequences() {
    
    if (userInput.inSequence.empty())
        return;

    std::vector<InSegment*> inSegments = *genome->getInSegments();
    
    for (InSegment *inSegment : inSegments) {
        lg.verbose("Processing segment: " + inSegment->getSeqHeader());
        DBGtoVariants(inSegment);
    }
}

bool DBG::DBGtoVariants(InSegment *inSegment) {
    
    Log threadLog;
    
    uint64_t explored = 0, len = inSegment->getSegmentLen();
    
    if (len<k)
        return true;
    
    uint64_t kcount = len-k+1;
    
    ParallelMap32* localGraphCache = new ParallelMap32;
    std::array<uint16_t, 2> mapRange;
    std::vector<std::deque<DBGpath>> variants;
    
    unsigned char* first = (unsigned char*)inSegment->getInSequencePtr()->c_str(); // prepare segment
    uint8_t* str = new uint8_t[len];
    bool* visited = new bool[len]{false};
    
    for (uint64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    while(explored < kcount) {

        mapRange = {0,0};
        std::deque<uint64_t> targetsQueue;
        phmap::parallel_flat_hash_map<uint64_t,bool> targetsMap;

        while (mapRange[1] < mapCount) {
            
            mapRange = computeMapRange(mapRange);
            loadMapRange(mapRange);
            
            uint64_t key, i;
            ParallelMap *map;
            ParallelMap32 *map32;
            bool isFw = false;
            
            for (uint16_t pos = 0; pos < userInput.maxSpan; ++pos) { // populate targets
                if (pos+k < kcount) {
                    key = hash(str+pos+k);
                    targetsQueue.push_back(key);
                    targetsMap[key];
                }
            }
            
            for (uint64_t c = 0; c<kcount; ++c){
                
                targetsMap.erase(targetsQueue.front()); // update targets
                targetsQueue.pop_front();
                if (c+k+userInput.maxSpan < kcount) {
                    key = hash(str+c+k+userInput.maxSpan);
                    targetsMap[key];
                    targetsQueue.push_back(key);
                }
                
                if(!visited[c]) {
                    
                    key = hash(str+c, &isFw);
                    i = key % mapCount;
                    
                    if (i >= mapRange[0] && i < mapRange[1]) {
                        
                        map = maps[i];
                        auto it = map->find(key);
                        std::pair<uint64_t,DBGkmer32> pair;
                        if (it != map->end()) {
                            
                            if (pair.second.cov == 255) {
                                map32 = maps32[i];
                                auto it = map32->find(key);
                                if (it == map32->end()) {
                                    std::cerr<<"Error: int32 map missing 255 value from int8 map"<<std::endl;
                                    exit(EXIT_FAILURE);
                                }
                                pair = *it;
                            }else{
                                pair = *it;
                            }
                            auto results = searchVariants(pair, isFw, hash(str+c+1, &isFw), mapRange, targetsQueue, targetsMap, localGraphCache);
                            explored += results.first;
                            if (results.first) {
                                for (DBGpath &path : results.second)
                                    path.pos = c+k;

                                if (results.second.size() != 0)
                                    variants.push_back(results.second);
                                
                                visited[c] = true;
                            }
                        }else{
                            explored += 1;
                            visited[c] = true;
                        }
                    }
                }
            }
            deleteMapRange(mapRange);
        }
    }
    delete localGraphCache;
    delete[] str;
    delete[] visited;
    
    inSegment->addVariants(variants);
        
    std::string ext = getFileExt("." + userInput.outFile);
    if (ext == "gfa" || ext == "gfa2" || ext == "gfa.gz" || ext == "gfa2.gz")
        variantsToGFA(inSegment, threadLog);
    
    std::unique_lock<std::mutex> lck(mtx);
    logs.push_back(threadLog);
    
    return true;
}

std::pair<bool,std::deque<DBGpath>> DBG::searchVariants(std::pair<const uint64_t,DBGkmer32> source, bool isSourceFw, uint64_t ref, std::array<uint16_t, 2> mapRange, const std::deque<uint64_t> &targetsQueue, const phmap::parallel_flat_hash_map<uint64_t,bool> &targetsMap, ParallelMap32* localGraphCache) { // dijkstra variant search
    
    bool explored = false; // true if we reached a node in the original graph
    std::vector<uint64_t> destinations;
    FibonacciHeap<std::pair<const uint64_t, DBGkmer32>*> Q; // node priority queue Q
    phmap::parallel_flat_hash_map<uint64_t,uint8_t> dist; // distance table
    phmap::parallel_flat_hash_map<uint64_t,std::pair<uint64_t,bool>> prev; // path table
    std::deque<DBGpath> discoveredPaths;
    
    dist[source.first] = 1;
    Q.insert(&source, 1); // associated priority equals dist[·]
    
    uint64_t key = source.first;
    int16_t depth = 0;
    bool direction = true, isFw;
    
    while (Q.size() > 0 && depth < userInput.kmerDepth + 1) { // the main loop
        explored = false; // if there are still node in the queue we cannot be done
        ParallelMap *map;
        //        ParallelMap32 *map32;
        
        std::pair<const uint64_t, DBGkmer32>* u = Q.extractMin(); // remove and return best vertex
        auto got = prev.find(u->first); // check direction
        if (got != prev.end()) {
            direction = got->second.second;
        }
        auto checkNext = [&,this] (uint64_t key, bool direction) {
            auto startNode = targetsMap.find(key);
            if (startNode == targetsMap.end()) { // if we connect to the original graph we are done
                auto nextKmer = localGraphCache->find(key); // check if the node is in the cache (already visited)
                
                if (nextKmer == localGraphCache->end()) { // we did not cache this node before
                    uint64_t m = key % mapCount;
                    if (m >= mapRange[0] && m < mapRange[1]) { // the node is in not cached but is available to visit now
                        map = maps[m];
                        auto got = map->find(key);
                        nextKmer = localGraphCache->insert(*got).first; // cache node for future iterations
                    }else{
                        return false;
                    }
                }
                uint8_t alt = dist[u->first]; // g(n)
                if (alt < std::numeric_limits<uint8_t>::max())
                    alt += 1; // Graph.Edges(u, v) << actual weight g(h)
                auto got = dist.find(nextKmer->first);
                if (got == dist.end()) { // if the next node was not seen before we add it to the queue
                    dist[nextKmer->first] = std::numeric_limits<uint8_t>::max(); // unknown distance from source to v
                    Q.insert(&*nextKmer, 0);
                }
                if (alt < dist[nextKmer->first]) {
                    prev[nextKmer->first] = std::make_pair(u->first,direction);
                    dist[nextKmer->first] = alt;
                    Q.decreaseKey(&*nextKmer, alt);
                }
            }
            return true;
        };
        uint8_t edgeCount = 0, exploredCount = 0;
        std::vector<std::tuple<uint64_t,bool,bool>> candidatePaths;
        
        for (uint8_t i = 0; i<4; ++i) { // first we collect all candidate paths
                
            if (depth == 0)
                direction = isSourceFw ? true : false;
            
            if (direction ? u->second.fw[i] : u->second.bw[i] > userInput.covCutOff) {
                uint8_t nextKmer[k];
                buildNextKmer(nextKmer, u->first, i, direction); // compute next node
                key = hash(nextKmer, &isFw);
                if (key != ref) { // we never rediscover the ref path
                    candidatePaths.push_back(std::make_tuple(key, isFw, direction));
                    ++edgeCount;
                }
            }
        }
        for (std::tuple<uint64_t,bool,bool> path : candidatePaths) { // for each path, including ref, we try to explore it
            uint64_t key = std::get<0>(path);
            bool isFw = std::get<1>(path);
            bool direction = std::get<2>(path);
            
            bool found = checkNext(key, isFw ? direction : !direction);
            if (found) { // if the next node exists
                ++exploredCount;
                if (targetsMap.find(key) != targetsMap.end()) { // if we reconnected to the ref path
                    prev[key] = std::make_pair(u->first,direction);
                    destinations.push_back(key); // we found a new alternate path
                }
            }
        }
        depth += 1;
        
        if(edgeCount == exploredCount || depth == userInput.kmerDepth + 1 || destinations.size() >= 10) // everything explored/found, depth reached, or top10
            explored = true;
    }
    if (destinations.size() > 0) { // traverse from target to source
        for (uint64_t destination : destinations) {
            DBGpath newPath;
            std::string endSequence = reverseHash(destination);
            uint16_t i = 0, refLen = std::find(targetsQueue.begin(), targetsQueue.end(), destination) - targetsQueue.begin() + k;
            uint64_t prevNode = prev[destination].first;

            while (prevNode != source.first) { // construct the shortest path with a stack S
                prevNode = prev[prevNode].first;
                ++i;
            }
            std::cout<<+i<<" "<<+refLen<<std::endl;
            prevNode = prev[destination].first;
            bool direction = prev[prevNode].second;
            int16_t b = i-refLen;
            if (refLen > k) {
                newPath.type = COM;
                newPath.refLen = refLen-k+1;
                b = refLen - k;
            }
            else if (i == refLen)
                newPath.type = SNV;
            else if (i > refLen) {
                newPath.type = DEL;
                --b;
                prevNode = prev[prevNode].first;
                direction = prev[prevNode].second;
            }
            else
                newPath.type = INS;

            while (b >= 0) {
                newPath.sequence.push_back(direction ? reverseHash(prevNode)[0] : revCom(reverseHash(prevNode)[k-1]));
                prevNode = prev[prevNode].first;
                direction = prev[prevNode].second;
                --b;
            }
            reverse(newPath.sequence.begin(), newPath.sequence.end());
            discoveredPaths.push_back(newPath);
        }
    }
    if (explored) { // we exahusted the search
        for (auto node : dist) // clear the cache for this source
            localGraphCache->erase(node.first);
    }
    return std::make_pair(explored,discoveredPaths);
}

bool DBG::variantsToGFA(InSegment *inSegment, Log &threadLog) {
    
    uint64_t processed = 0;
    uint32_t segmentCounter = 0, edgeCounter = 0, sUId, sUIdNew;
    std::string *oldSequence = inSegment->getInSequencePtr();
    std::string sHeader = inSegment->getSeqHeader();
    std::string* inSequence;
    Sequence* newSequence;
    InSegment* newSegment;
    std::vector<uint32_t> sUIds;
    uint32_t seqPos = inSegment->getSeqPos(), sId = 0, eId = 0;
    std::vector<std::deque<DBGpath>>& variants = inSegment->getVariants();
    
    for (std::deque<DBGpath> DBGpaths : variants) {
        
        threadLog.add("Introducing variants at pos: " + std::to_string(DBGpaths[0].pos+1));
        
        inSequence = new std::string(oldSequence->substr(processed, DBGpaths[0].pos-processed));
        newSequence = new Sequence{sHeader + "." + std::to_string(++segmentCounter), "", inSequence, NULL, seqPos};
        threadLog.add("Previous sequence: " + newSequence->header);
        newSegment = genome->traverseInSegment(newSequence, std::vector<Tag>(), sId++);
        sUId = newSegment->getuId();
        
        for (uint32_t sUIdprev : sUIds) {
            InEdge edge(seqPos, eId++, sUIdprev, sUId, '+', '+', "0M", sHeader + ".edge." + std::to_string(++edgeCounter));
            genome->appendEdge(edge);
        }
        sUIds.clear();
        
        uint32_t altCounter = 0;
        bool originalAdded = false;
        processed = DBGpaths[0].pos;
        
        for (DBGpath variant : DBGpaths) {
            
            if (variant.type != DEL && !originalAdded) {
                
                inSequence = new std::string(oldSequence->substr(DBGpaths[0].pos, 1));
                newSequence = new Sequence{sHeader + "." + std::to_string(++segmentCounter), "Candidate sequence", inSequence, NULL, seqPos};
                threadLog.add("Old segment from SNV/DEL error: " + newSequence->header + "\t" + *inSequence);
                newSegment = genome->traverseInSegment(newSequence, std::vector<Tag>(), sId++);
                sUIdNew = newSegment->getuId();
                sUIds.push_back(sUIdNew);
                
                InEdge edge(seqPos, eId++, sUId, sUIdNew, '+', '+', "0M", sHeader + ".edge." + std::to_string(++edgeCounter));
                genome->appendEdge(edge);
                
                originalAdded = true;
                ++processed;
                
            }
            if (variant.type == SNV || variant.type == DEL) {
                
                inSequence = new std::string(variant.sequence);
                newSequence = new Sequence{sHeader + "." + std::to_string(segmentCounter) + ".alt" + std::to_string(++altCounter), "Candidate sequence", inSequence, NULL, seqPos};
                threadLog.add("New segment from SNV error: " + newSequence->header + "\t" + *inSequence);
                newSegment = genome->traverseInSegment(newSequence, std::vector<Tag>(), sId++);
                sUIdNew = newSegment->getuId();
                sUIds.push_back(sUIdNew);
            }
            if (variant.type == SNV) {

                InEdge edge(seqPos, eId++, sUId, sUIdNew, '+', '+', "0M", sHeader + ".edge." + std::to_string(++edgeCounter));
                genome->appendEdge(edge);
                
            }else if (variant.type == INS) {
                
                sUIds.push_back(sUId);
                
            }else if (variant.type == DEL) {
                
                InEdge edge(seqPos, eId++, sUId, sUIdNew, '+', '+', "0M", sHeader + ".edge." + std::to_string(++edgeCounter));
                genome->appendEdge(edge);
                sUIds.push_back(sUId);
            }
        }
    }
    if (variants.size() > 0) { // residual sequence
        
        inSequence = new std::string(oldSequence->substr(processed));
        newSequence = new Sequence{sHeader + "." + std::to_string(++segmentCounter), "", inSequence, NULL, seqPos};
        threadLog.add("Previous sequence: " + newSequence->header);
        newSegment = genome->traverseInSegment(newSequence, std::vector<Tag>(), sId++);
        sUId = newSegment->getuId();
        
        for (uint32_t sUIdprev : sUIds) {
            InEdge edge(seqPos, eId++, sUIdprev, sUId, '+', '+', "0M", sHeader + ".edge." + std::to_string(++edgeCounter));
            genome->appendEdge(edge);
        }
        genome->deleteSegment(sHeader);
    }
    return true;
}

bool DBG::detectAnomalies(InSegment *inSegment, std::vector<uint64_t> &anomalies) { // find disconnected kmers in segments
    
    Log threadLog;
    threadLog.setId(inSegment->getuId());
        
    std::string sHeader = inSegment->getSeqHeader();
    ParallelMap *map;
    uint64_t key, i;
    bool isFw = false, anomaly = false;
        
    uint64_t len = inSegment->getSegmentLen();
    
    if (len<k)
        return true;
    
    uint64_t kcount = len-k+1;
    
    unsigned char* first = (unsigned char*)inSegment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];
    
    for (uint64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    for (uint64_t c = 0; c<kcount; ++c){

        key = hash(str+c, &isFw);
        i = key % mapCount;

        map = maps[i];
        auto got = map->find(key);
            
        // check for DBG consistency
        if (got != map->end()) {
            DBGkmer &dbgkmer = got->second;
            if (c < kcount-1 && ((isFw && dbgkmer.fw[*(str+c+k)] == 0) || (!isFw && dbgkmer.bw[3-*(str+c+k)] == 0))) // find alternative paths
                anomaly = true;
        }else{anomaly = true;}
        if (anomaly) {
            anomalies.push_back(c+k);
            threadLog.add("Anomaly at:\t" + sHeader + "\t" + std::to_string(c+k+1));
            anomaly = false;
        }
    }
    delete[] str;
    
    std::unique_lock<std::mutex> lck(mtx);
    logs.push_back(threadLog);
    
    return true;
    
}
