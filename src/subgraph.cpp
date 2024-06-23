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

#include "subgraph.h"

// help functions

void DBG::mergeSubgraphs() {
    
    for (ParallelMap32color *map1 : DBGTmpSubgraphs) {
        unionSum(map1, DBGsubgraph);
        delete map1;
    }
}

template<typename MAPTYPE>
bool DBG::mergeSubMaps(MAPTYPE* map1, MAPTYPE* map2, uint8_t subMapIndex) {
    
    auto& inner = map1->get_inner(subMapIndex);   // to retrieve the submap at given index
    auto& submap1 = inner.set_;        // can be a set or a map, depending on the type of map1
    auto& inner2 = map2->get_inner(subMapIndex);
    auto& submap2 = inner2.set_;
    
    for (auto pair : submap1) { // for each element in map1, find it in map2 and increase its value
        
        auto got = submap2.find(pair.first); // insert or find this kmer in the hash table
        if (got == submap2.end()){
            submap2.insert(pair);
        }else{

            DBGkmer32& dbgkmerMap = got->second;
        
            for (uint64_t w = 0; w<4; ++w) { // update weights
                
                if (LARGEST - dbgkmerMap.fw[w] >= pair.second.fw[w])
                    dbgkmerMap.fw[w] += pair.second.fw[w];
                else
                    dbgkmerMap.fw[w] = LARGEST;
                if (LARGEST - dbgkmerMap.bw[w] >= pair.second.bw[w])
                    dbgkmerMap.bw[w] += pair.second.bw[w];
                else
                    dbgkmerMap.bw[w] = LARGEST;
                
            }
            
            if (LARGEST - dbgkmerMap.cov >= pair.second.cov)
                dbgkmerMap.cov += pair.second.cov; // increase kmer coverage
            else
                dbgkmerMap.cov = LARGEST;
            
        };
        
    }
    
    return true;
    
}

template<typename MAPTYPE>
bool DBG::unionSum(MAPTYPE* map1, MAPTYPE* map2) {
    
    std::vector<std::function<bool()>> jobs;
    
    if (map1->subcnt() != map2->subcnt()) {
        fprintf(stderr, "Maps don't have the same numbers of submaps (%zu != %zu). Terminating.\n", map1->subcnt(), map2->subcnt());
        exit(EXIT_FAILURE);
    }
    
    for(std::size_t subMapIndex = 0; subMapIndex < map1->subcnt(); ++subMapIndex)
        jobs.push_back([this, map1, map2, subMapIndex] { return this->mergeSubMaps(map1, map2, subMapIndex); });
    
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
    return true;
    
}

// subgraph functions

void DBG::subgraph() {
    
    if (userInput.inSequence.empty())
        return;
    
    lg.verbose("Subsetting graph");
    
    if (!userInput.inBedInclude.empty() || userInput.pipeType == 'i') {
        
        StreamObj streamObj;
        std::shared_ptr<std::istream> stream = streamObj.openStream(userInput, 'i');
        std::string line, bedHeader;
        
        while (getline(*stream, line)) {
            
            uint64_t begin = 0, end = 0;
            std::istringstream iss(line);
            iss >> bedHeader >> begin >> end;
            userInput.bedIncludeList.pushCoordinates(bedHeader, begin, end);
        }
        lg.verbose("Finished reading BED include list");
        userInput.bedIncludeList = BEDPathsToSegments();
    }
    
    std::vector<std::function<bool()>> jobs;
    std::array<uint16_t, 2> mapRange = {0,0};

    while (mapRange[1] < mapCount) {

        mapRange = computeMapRange(mapRange);
        loadMapRange(mapRange);
        
        std::vector<InSegment*> inSegments = *genome->getInSegments();
        for (InSegment *inSegment : inSegments)
            jobs.push_back([this, inSegment, mapRange] { return DBGsubgraphFromSegment(inSegment, mapRange); });
        
        threadPool.queueJobs(jobs);
        jobWait(threadPool);
        jobs.clear();
        
        deleteMapRange(mapRange);

    }
    lg.verbose("Merging subgraphs");
    mergeSubgraphs();
    lg.verbose("Searching graph");
    if (userInput.travAlgorithm == "best-first") {
        if (userInput.kmerDepth == -1)
            userInput.kmerDepth = userInput.kmerLen; // unidirectional search
        bestFirst();
    }else if (userInput.travAlgorithm == "traversal") {
        if (userInput.kmerDepth == -1)
            userInput.kmerDepth = std::ceil((float)userInput.kmerLen/2); // kmer search is in both directions
        traversal();
    }else{
        fprintf(stderr, "Cannot find input algorithm (%s). Terminating.\n", userInput.travAlgorithm.c_str());
        exit(EXIT_FAILURE);
    }
    lg.verbose("Computing summary graph");
    summary(*DBGsubgraph);
    lg.verbose("Generating GFA");
    DBGgraphToGFA();
    
}

void DBG::summary(ParallelMap32color& DBGsubgraph) {
    
    uint64_t tot = 0, kmersUnique = 0, kmersDistinct = DBGsubgraph.size(), edgeCount = 0;
    phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;

    for (auto pair : DBGsubgraph) {
        
        if (pair.second.cov == 1)
            ++kmersUnique;
        
        for (uint8_t w = 0; w<4; ++w) // update weights
            edgeCount += pair.second.fw[w] > 0 ? 1 : 0 + pair.second.bw[w] > 0 ? 1 : 0;
        
        ++hist[pair.second.cov];
    }
    for (auto pair : hist)
        tot += pair.first * pair.second;
        
    uint64_t missing = pow(4,k)-kmersDistinct;
    std::cout<<"Subgraph summary statistics:\n"
             <<"Total kmers: "<<tot<<"\n"
             <<"Unique kmers: "<<kmersUnique<<"\n"
             <<"Distinct kmers: "<<kmersDistinct<<"\n"
             <<"Missing kmers: "<<missing<<"\n"
             <<"Total edges: "<<edgeCount<<"\n";
}

bool DBG::DBGsubgraphFromSegment(InSegment *inSegment, std::array<uint16_t, 2> mapRange) {
    
    Log threadLog;
    threadLog.setId(inSegment->getuId());
        
    std::string sHeader = inSegment->getSeqHeader();
    ParallelMap *map;
    ParallelMap32 *map32;
    ParallelMap32color *segmentSubmap = new ParallelMap32color;
    uint64_t key, i;
    bool isFw = false;
    std::vector<uint64_t> segmentCoordinates;
    
    uint64_t len = inSegment->getSegmentLen();
    if (len<k)
        return true;

    uint64_t kcount = len-k+1;
    
    unsigned char* first = (unsigned char*)inSegment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];
    
    for (uint64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    std::vector<std::pair<uint64_t,uint64_t>> span = {std::make_pair(0, kcount)};
    
    if (userInput.inBedInclude != "") {
        
        auto coordinates = userInput.bedIncludeList.getCoordinates();
        auto got = coordinates.find(inSegment->getSeqHeader());
        if (got != coordinates.end())
            span = got->second;
        else
            span.clear();
    }
    
    for (auto coordinate : span) {
        for (uint64_t p = coordinate.first; p < coordinate.second; ++p) {
            
            key = hash(str+p, &isFw);
            i = key % mapCount;
            
            if (i >= mapRange[0] && i < mapRange[1]) {
                
                map = maps[i];
                auto got = map->find(key);
                
                if (got != map->end()) {
                    if (got->second.cov != 255) {
                        DBGkmer32color dbgKmer32color(got->second);
                        dbgKmer32color.color = 1;
                        segmentSubmap->insert(std::make_pair(got->first,dbgKmer32color));
                    }else{
                        map32 = maps32[i];
                        auto got = map32->find(key);
                        DBGkmer32color dbgKmer32color(got->second);
                        dbgKmer32color.color = 1;
                        segmentSubmap->insert(std::make_pair(got->first,dbgKmer32color));
                    }
                }else if(!userInput.noReference){ // construct the kmer
                    
                    DBGkmer32color dbgKmer32color;
                    dbgKmer32color.color = 2;
                    edgeBit edges;
                    
                    if (isFw){
                        if (ctoi[*(first+p+k)] <= 3)
                            edges.assign(ctoi[*(first+p+k)]);
                        if (p > 0 && *(str+p-1) <= 3)
                            edges.assign(4+*(str+p-1));
                    }else{
                        if (p > 0 && *(str+p-1) <= 3)
                            edges.assign(3-*(str+p-1));
                        if (ctoi[*(first+p+k)] <= 3)
                            edges.assign(4+3-ctoi[*(first+p+k)]);
                    }
                    
                    for (uint64_t w = 0; w<4; ++w) { // update weights
                       
                        dbgKmer32color.fw[w] += edges.read(w);
                        dbgKmer32color.bw[w] += edges.read(4+w);
                    }
                    if (dbgKmer32color.cov < LARGEST)
                        ++dbgKmer32color.cov; // increase kmer coverage
                    
                    segmentSubmap->insert(std::make_pair(key,dbgKmer32color));
                }
            }
        }
    }
    delete[] str;
    
    std::unique_lock<std::mutex> lck (mtx);
    DBGTmpSubgraphs.push_back(segmentSubmap);
    logs.push_back(threadLog);
    
    return true;
}

void DBG::traversal() {
    
    ParallelMap32color candidates, newCandidates;
    ParallelMap32color* subgraph = DBGsubgraph;
    
    std::array<uint16_t, 2> mapRange = {0,0};
    for (uint8_t i = 0; i < userInput.kmerDepth; ++i) {

        mapRange = {0,0};

        while (mapRange[1] < mapCount) {

            mapRange = computeMapRange(mapRange);
            loadMapRange(mapRange);
            newCandidates = traversalPass(subgraph, mapRange);
            deleteMapRange(mapRange);
            candidates.insert(newCandidates.begin(), newCandidates.end());
            subgraph = &newCandidates;
        }
    }
    DBGsubgraph->insert(candidates.begin(), candidates.end());
}

ParallelMap32color DBG::traversalPass(ParallelMap32color* subgraph, std::array<uint16_t, 2> mapRange) {
    
    ParallelMap32color newCandidates;
    
    for (auto pair : *subgraph) {
        
        for (uint8_t i = 0; i<4; ++i) { // forward edges
            if (pair.second.fw[i] != 0) {
                
                uint8_t nextKmer[k];
                std::string firstKmer = reverseHash(pair.first);
                firstKmer.push_back(itoc[i]);
                for (uint8_t e = 0; e<k; ++e)
                    nextKmer[e] = ctoi[(unsigned char)firstKmer[e+1]];
                
                ParallelMap *map;
                ParallelMap32 *map32;
                bool isFw = false;
                uint64_t key = hash(nextKmer, &isFw);
                
                auto got = DBGsubgraph->find(key); // avoid trying implicit searches
                
                if (got == DBGsubgraph->end()) {
                    
                    uint64_t m = key % mapCount;
                    
                    if (m >= mapRange[0] && m < mapRange[1]) {
                        
                        map = maps[m];
                        auto got = map->find(key);
                        
                        if (got != map->end()) {
                            
                            if (got->second.cov != 255) {
                                DBGkmer32color dbgKmer32color(got->second);
                                newCandidates.insert(std::make_pair(got->first,dbgKmer32color));
                            }else{
                                map32 = maps32[m];
                                auto got = map32->find(key);
                                DBGkmer32color dbgKmer32color(got->second);
                                newCandidates.insert(std::make_pair(got->first,dbgKmer32color));
                            }
                        }
                    }
                }
            }
        }
        for (uint8_t i = 0; i<4; ++i) { // reverse edges
            if (pair.second.bw[i] != 0) {
                
                uint8_t nextKmer[k];
                std::string firstKmer;
                firstKmer.push_back(itoc[i]);
                firstKmer.append(reverseHash(pair.first));
                
                for (uint8_t e = 0; e<k; ++e)
                    nextKmer[e] = ctoi[(unsigned char)firstKmer[e]];
                
                ParallelMap *map;
                ParallelMap32 *map32;
                bool isFw = false;
                uint64_t key = hash(nextKmer, &isFw);
                
                auto got = DBGsubgraph->find(key); // avoid trying implicit searches
                
                if (got == DBGsubgraph->end()) {
                    
                    uint64_t m = key % mapCount;
                    
                    if (m >= mapRange[0] && m < mapRange[1]) {
                        
                        map = maps[m];
                        auto got = map->find(key);
                        
                        if (got != map->end()) {
                            if (got->second.cov != 255) {
                                DBGkmer32color dbgKmer32color(got->second);
                                newCandidates.insert(std::make_pair(got->first,dbgKmer32color));
                            }else{
                                map32 = maps32[m];
                                auto got = map32->find(key);
                                DBGkmer32color dbgKmer32color(got->second);
                                newCandidates.insert(std::make_pair(got->first,dbgKmer32color));
                            }
                        }
                    }
                }
            }
        }
    }
    return newCandidates;
}

void DBG::bestFirst() {
    
    ParallelMap32color* candidates = new ParallelMap32color;
    uint32_t explored = 0, total = DBGsubgraph->size();
    std::array<uint16_t, 2> mapRange;
    ParallelMap32color* DBGsubgraphCpy = new ParallelMap32color;
    while(explored < total) {

        mapRange = {0,0};

        while (mapRange[1] < mapCount) {

            mapRange = computeMapRange(mapRange);
            loadMapRange(mapRange);
            for (auto pair : *DBGsubgraph) {
                auto results = dijkstra(pair, mapRange);;
                explored += results.first;
                if (results.first) {
                    candidates->insert(results.second.begin(), results.second.end());
                    DBGsubgraphCpy->insert(pair);
//                    DBGsubgraph->erase(pair.first);
                }
//                std::cout<<DBGsubgraphCpy.size()<<std::endl;
            }
            deleteMapRange(mapRange);
        }
    }
    DBGsubgraphCpy->insert(candidates->begin(), candidates->end());
    delete DBGsubgraph;
    DBGsubgraph = DBGsubgraphCpy;
    delete candidates;
}

std::pair<bool,ParallelMap32color> DBG::dijkstra(std::pair<uint64_t,DBGkmer32color> source, std::array<uint16_t, 2> mapRange) {
    
    bool explored = false; // true if we reached a node in the original graph
    std::vector<uint64_t> destinations;
    FibonacciHeap<std::pair<const uint64_t, DBGkmer32>*> Q; // node priority queue Q
    phmap::parallel_flat_hash_map<uint64_t,uint8_t> dist;
    phmap::parallel_flat_hash_map<uint64_t,uint64_t> prev; // distance table
    ParallelMap32color discoveredNodes;
    
    dist[source.first] = 1;
    std::pair<const uint64_t, DBGkmer32> firstKmer = std::make_pair(source.first,source.second);
    Q.insert(&firstKmer, 1); // associated priority equals dist[·]
    
    uint64_t key;
    int16_t depth = 0;
    
    while (Q.size() > 0 && depth < userInput.kmerDepth + 1) { // The main loop
        ParallelMap *map;
        //        ParallelMap32 *map32;
        
        bool isFw = false;
        std::pair<const uint64_t, DBGkmer32>* u = Q.extractMin(); // Remove and return best vertex

        auto checkNext = [&,this] (uint64_t key) {
            auto startNode = DBGsubgraph->find(key);
            if (startNode == DBGsubgraph->end()) { // if we connect to the original graph we are done
                auto nextKmer = graphCache->find(key); // check if the node is in the cache (already visited)
                
                if (nextKmer == graphCache->end()) { // we cached this node before
                    uint64_t m = key % mapCount;
                    if (m >= mapRange[0] && m < mapRange[1]) { // the node is in not cached but is available to visit now
                        map = maps[m];
                        auto got = map->find(key);
                        nextKmer = graphCache->insert(*got).first; // cache node for future iterations
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
                    prev[nextKmer->first] = u->first;
                    dist[nextKmer->first] = alt;
                    Q.decreaseKey(&*nextKmer, alt);
                }
            }
            return true;
        };
        uint8_t edgeCount = 0, exploredCount = 0;
        for (uint8_t i = 0; i<4; ++i) { // forward edges
            if (u->second.fw[i] > userInput.covCutOff) {
                uint8_t nextKmer[k];
                buildNextKmer(nextKmer, u->first, i, true); // compute next node
                key = hash(nextKmer, &isFw);
                bool found = checkNext(key);
                if (found) {
                    ++exploredCount;
                    if (depth > 1 && DBGsubgraph->find(key) != DBGsubgraph->end())
                        destinations.push_back(u->first);
                }
                ++edgeCount;
            }
            if (u->second.bw[i] > userInput.covCutOff) { // backward edges
                uint8_t nextKmer[k];
                buildNextKmer(nextKmer, u->first, i, false); // compute next node
                key = hash(nextKmer, &isFw);
                bool found = checkNext(key);
                if (found) {
                    ++exploredCount;
                    if (depth > 1 && DBGsubgraph->find(key) != DBGsubgraph->end())
                        destinations.push_back(u->first);
                }
                ++edgeCount;
            }
        }
        depth += 1;
        if(edgeCount == exploredCount || depth == userInput.kmerDepth + 1 || destinations.size() >= 10) // everything explored/found, depth reached, or top10
            explored = true;
    }
    if (destinations.size() > 0) { // traverse from target to source
        for (uint64_t destination : destinations) {
            while (destination != source.first) { // construct the shortest path with a stack S
                discoveredNodes.insert(*graphCache->find(destination)); // push the vertex onto the stack
                dist.erase(destination);
                destination = prev[destination];
            }
        }
    }else{
        
    }
    if (explored) { // we exahusted the search
        for (auto node : dist) // clear the cache for this source
            graphCache->erase(node.first);
    }
    return std::make_pair(explored,discoveredNodes);
}

void DBG::buildNextKmer(uint8_t* nextKmer, uint64_t hash, uint8_t nextBase, bool fw) {
    
    std::string firstKmer;
    
    if (fw) {
        firstKmer.append(reverseHash(hash));
        firstKmer.push_back(itoc[nextBase]);
        for (uint8_t e = 0; e<k; ++e)
            nextKmer[e] = ctoi[(unsigned char)firstKmer[e+1]];
    }else{
        firstKmer.push_back(itoc[nextBase]);
        firstKmer.append(reverseHash(hash));
        for (uint8_t e = 0; e<k; ++e)
            nextKmer[e] = ctoi[(unsigned char)firstKmer[e]];
        
    }
}










