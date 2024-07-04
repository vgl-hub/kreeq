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

double errorRate(uint64_t missingKmers, uint64_t totalKmers, uint8_t k){ // estimate QV from missing kmers
    
    return 1 - pow(1 - (double) missingKmers/totalKmers, (double) 1/k);
    
}

void DBG::loadGenome(InSequencesDBG *genome) {
    
    this->genome = genome;
}

void DBG::validateSequences() {
    
    if (userInput.inSequence.empty())
        return;
    
    lg.verbose("Validating sequence");
    
    std::vector<std::function<bool()>> jobs;
    std::vector<InSegment*> *segments = genome->getInSegments();
    std::array<uint16_t, 2> mapRange = {0,0};
    genome->generateValidationVector();
    
    while (mapRange[1] < mapCount) {
        
        mapRange = computeMapRange(mapRange);
        
        loadMapRange(mapRange);
        
        for (uint32_t s = 0; s < segments->size(); ++s)
            jobs.push_back([this, s, mapRange] { return evaluateSegment(s, mapRange); });
        
        threadPool.queueJobs(jobs);
        
        jobWait(threadPool);
        
        jobs.clear();
            
        deleteMapRange(mapRange);
        
    }
    
    if (userInput.outFile.find(".") != std::string::npos || userInput.outFile == "") {
        
        std::cout<<"Missing"<<"\t"
        <<"Total"<<"\t"
        <<"QV"<<"\t"
        <<"Error"<<"\t"
        <<"k"<<"\t"
        <<"Method"<<std::endl;
        
        double merquryError = errorRate(totMissingKmers, totKcount, k), merquryQV = -10*log10(merquryError);
        
        std::cout<<totMissingKmers<<"\t"
        <<totKcount<<"\t"
        <<merquryQV<<"\t"
        <<merquryError<<"\t"
        <<std::to_string(k)<<"\t"
        <<"Merqury"<<std::endl;
        
        double kreeqError = errorRate(totMissingKmers+totEdgeMissingKmers, totKcount, k), kreeqQV = -10*log10(kreeqError);
        
        std::cout
        <<totMissingKmers+totEdgeMissingKmers<<"\t"
        <<totKcount<<"\t"
        <<kreeqQV<<"\t"
        <<kreeqError<<"\t"
        <<std::to_string(k)<<"\t"
        <<"Kreeq"<<std::endl;
        
    }
    
}

bool DBG::evaluateSegment(uint32_t s, std::array<uint16_t, 2> mapRange) {
    
    std::vector<uint64_t> missingKmers;
    std::vector<uint64_t> edgeMissingKmers;
    
    std::vector<InSegment*> *segments = genome->getInSegments();
    std::vector<DBGbase*> *dbgbases = genome->getInSegmentsDBG();
    
    InSegment *segment = (*segments)[s];
    DBGbase *DBGsequence = (*dbgbases)[s];
    
    uint64_t len = segment->getSegmentLen();
    
    if (len<k)
        return true;
    
    uint64_t kcount = len-k+1, kmers = 0;
    
    unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];
    
    for (uint64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    uint64_t key, i;
    
    ParallelMap *map;
    ParallelMap32 *map32;
    
    // kreeq QV
    bool isFw = false;
    
    // std::cout<<segment->getSeqHeader()<<std::endl;
    for (uint64_t c = 0; c<kcount; ++c){
        
        key = hash(str+c, &isFw);
        i = key % mapCount;
        
//        std::cout<<"\n"<<itoc[*(str+c)]<<"\t"<<c<<"\t"<<isFw<<std::endl;
        
        if (i >= mapRange[0] && i < mapRange[1]) {
            
            map = maps[i];
            auto it = map->find(key);
            DBGkmer32 khmer;
            if (it != map->end()) {
                khmer = it->second;
                
                if (khmer.cov == 255) {
                    map32 = maps32[i];
                    auto it = map32->find(key);
                    if (it == map32->end()) {
                        std::cerr<<"Error: int32 map missing 255 value from int8 map"<<std::endl;
                        exit(EXIT_FAILURE);
                    }
                    khmer = it->second;
                }
                
                DBGsequence[c].cov = khmer.cov;
                DBGsequence[c].isFw = isFw;
            }
            
            if (DBGsequence[c].cov == 0) // merqury QV
                missingKmers.push_back(c);
            else if (DBGsequence[c].cov < userInput.covCutOff) // merqury QV with cutoff
                missingKmers.push_back(c);
            else { // kreeq QV
                bool noEdgeLeft = false, noEdgeRight = false;
                if (DBGsequence[c].isFw){

                    if (c<kcount-1) {
                        if (khmer.fw[*(str+c+k)] != 0)
                            DBGsequence[c].fw = khmer.fw[*(str+c+k)];
                        else
                            noEdgeRight = true;
                    }
                    if (c>0) {
                        if (khmer.bw[*(str+c-1)] != 0)
                            DBGsequence[c].bw = khmer.bw[*(str+c-1)];
                        else
                            noEdgeLeft = true;
                        // std::cout<<"edge error1"<<std::endl;
                    }
                    
                }else{
                    
                    if (c>0) {
                        if (khmer.fw[3-*(str+c-1)] != 0)
                            DBGsequence[c].fw = khmer.fw[3-*(str+c-1)];
                        else
                            noEdgeLeft = true;
                    }
                    
                    if (c<kcount-1) {
                        if (khmer.bw[3-*(str+c+k)] != 0)
                            DBGsequence[c].bw = khmer.bw[3-*(str+c+k)];
                        else
                            noEdgeRight = true;
                        // std::cout<<"edge error2"<<std::endl;
                    }
                }
                if (noEdgeLeft && noEdgeRight) {
                    edgeMissingKmers.push_back(c);
//                    std::cout<<"edge error"<<std::endl;
                }
            }
            ++kmers;

        }
    }
    
    delete[] str;
    
    totMissingKmers += missingKmers.size();
    totKcount += kmers;
    totEdgeMissingKmers += edgeMissingKmers.size();
    
    return true;
    
}

std::pair<DBGkmer*,bool> DBG::findDBGkmer(uint8_t *origin) {
    
    uint64_t key, i;
    bool isFw = false;
    key = hash(origin, &isFw);
    i = key % mapCount;
    ParallelMap *map = maps[i];
    auto got = map->find(key);
    if (got != map->end())
        return std::make_pair(&(got->second), isFw);
    else
        return std::make_pair(nullptr,false);
    
}

int8_t DBG::checkNext(uint8_t *currentKmer, uint8_t *nextBase) {
    
    int8_t score = 0;
    uint8_t nextKmer[k];
    memcpy(nextKmer, currentKmer, k);
    
    for (uint8_t i = 0; i < k; ++i) {
        
        memmove(nextKmer, nextKmer+1, k-1 * sizeof(uint8_t));
        nextKmer[k-1] = *(nextBase+i);
            
        auto it = findDBGkmer(nextKmer);
        DBGkmer *dbgOrigin = it.first;
        
        if (dbgOrigin == nullptr)
            --score;
        else
            ++score;
        
    }
    
    return score;
    
}

bool DBG::loadSegmentCoordinates(InSegment *inSegment, std::vector<uint64_t> &segmentCoordinates) {
    
    auto coordinates = userInput.bedIncludeList.getCoordinates();
    auto got = coordinates.find(inSegment->getSeqHeader());
    if (got != coordinates.end()) {
        for (auto coordinate : got->second) {
            for (uint64_t i = coordinate.first; i < coordinate.second; ++i)
                segmentCoordinates.push_back(i);
        }
    }
    return true;
}

BedCoordinates DBG::BEDPathsToSegments() { // project path coordinates to segments
    
    BedCoordinates inBedIncludeSegments;
    
    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    auto coordinates = userInput.bedIncludeList.getCoordinates();
    
    for (InPath& path : inPaths) {
        
        auto got = coordinates.find(path.getHeader());
        if (got == coordinates.end())
            continue;
        
        unsigned int cUId = 0, gapLen = 0;
        std::vector<PathComponent> pathComponents = path.getComponents();
        uint64_t absPos = 0;
        std::vector<std::pair<uint64_t,uint64_t>>::iterator pathCoordinatePtr = got->second.begin();
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            
            if (component->componentType == SEGMENT) {
                
                auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                if (component->orientation == '+') {
                    
                    while (pathCoordinatePtr->first > absPos && pathCoordinatePtr->first < absPos + (*inSegment)->getSegmentLen()) {
                        uint64_t cBegin = pathCoordinatePtr->first - absPos, cEnd = pathCoordinatePtr->second - absPos;
                        inBedIncludeSegments.pushCoordinates((*inSegment)->getSeqHeader(), cBegin, cEnd);
                        ++pathCoordinatePtr;
                        if (pathCoordinatePtr == got->second.end())
                            break;
                    }
                }else{} // GFA not handled yet
            }else if (component->componentType == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                gapLen = inGap->getDist(component->start - component->end);
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
        }
    }
    
    return inBedIncludeSegments;
}



std::string colorPalette(uint8_t value){
    
    const static phmap::parallel_flat_hash_map<uint8_t,std::string> value_to_color{ // map each value to a color
        {0,"gray"},
        {1,"blue"},
        {2,"red"},
    };
    if (value > value_to_color.size()) {
        fprintf(stderr, "Node color value %i does not exist\n", value);
        exit(EXIT_FAILURE);
    }
    return (value_to_color.find(value)->second);
}

void DBG::nextKmerFromString(uint8_t *nextKmer, std::string *sequence, uint64_t start, uint8_t nextBase){
    
    std::string kmer = sequence->substr(start, k);
    kmer.push_back(itoc[nextBase]); // append its base
    for (uint8_t e = 0; e<k; ++e)
        nextKmer[e] = ctoi[(unsigned char)kmer[e+1]];
    
}

void DBG::collapseNodes() {
    
    uint32_t idCounter = 0, seqPos = 0, edgeCounter = 0;
    phmap::flat_hash_map<std::string, unsigned int>& headersToIds = *GFAsubgraph.getHash1();
    phmap::parallel_flat_hash_map<uint64_t, std::tuple<DBGkmer32,uint32_t,bool>> residualEdges; // hash, kmer, G' node, node side
    
    auto extend = [&,this] (std::string &seed, int8_t direction) {
        
        uint64_t key, baseCounter = 0;
        bool isFw;
        uint8_t nextKmer[k];
        for (uint8_t e = 0; e<k; ++e)
            nextKmer[e] = ctoi[(unsigned char)seed[e]];
        key = hash(nextKmer, &isFw);
        
        std::pair<uint64_t, DBGkmer32color> node = *DBGsubgraph->find(key);
        
        if ((direction ? node.second.fwCount() : node.second.bwCount()) > 1) {
            std::cout<<"Branching node side, cannot extend. Terminating."<<std::endl;
            exit(EXIT_FAILURE);
        }else if ((direction ? node.second.fwCount() : node.second.bwCount()) == 0){
            std::cout<<"Dead end, cannot extend. Terminating."<<std::endl;
            exit(EXIT_FAILURE);
        }
        
        while (true) {
            
            uint8_t i = isFw ? node.second.fwEdgeIndexes()[0] : 3-node.second.bwEdgeIndexes()[0];
            
            uint8_t nextKmer[k];
            nextKmerFromString(nextKmer, &seed, baseCounter++, i);
            key = hash(nextKmer, &isFw);
            
//                std::cout<<seed<<std::endl;
//                std::cout<<+i<<std::endl;
            auto prevNode = node;
            auto got = DBGsubgraph->find(key);
            
            if (got == DBGsubgraph->end()) { // we found a novel dead end
                auto got = residualEdges.find(key); // we couldn't find the node as it was already visited and deleted
                if(got != residualEdges.end()) {
                    residualEdges[prevNode.first] = std::make_tuple(prevNode.second,idCounter,direction);
                    
                }
                break; // real dead ends
            }

            node = *got;

            std::vector<uint8_t> nextFrontEdges = isFw ? node.second.fwEdgeIndexes() : node.second.bwEdgeIndexes();
            std::vector<uint8_t> nextBackEdges = isFw ? node.second.bwEdgeIndexes() : node.second.fwEdgeIndexes();
            
            if(nextBackEdges.size() > 1)  { // this node is branching back, we cannot include it
                residualEdges[prevNode.first] = std::make_tuple(prevNode.second,idCounter,direction);
                break;
            }
            
            seed.push_back(itoc[i]); // append its base
            DBGsubgraph->erase(key); // we can now safely erase as these nodes don't need to be stored
            
            if (nextFrontEdges.size() == 0) { // we found a dead end
                break;
            } else if (nextFrontEdges.size() > 1) { // we found a potentially fw branching node, if true, nothing more to be done
                residualEdges[node.first] = std::make_tuple(node.second,idCounter,direction); // we preserve the edge information
                break;
            }
        }
    };
    
    while (DBGsubgraph->size() != 0) { // until all nodes have been merged or outputted
        
        auto pair = DBGsubgraph->begin(); // pick a random node
        std::string frontSequence = reverseHash(pair->first); // we grow the sequence in both directions
        std::string backSequence = revCom(reverseHash(pair->first));
    
        uint8_t edgeCounts[2] = {pair->second.bwCount(), pair->second.fwCount()};
        
        if (edgeCounts[0] == 1 || edgeCounts[1] == 1) { // we are in the middle or at a partial branch, we can start merging
                
            for (int8_t direction = 1; direction >= 0; --direction) { // we potentially extend in both directions
                
                if (edgeCounts[direction] == 1) { // we can extend if we are at a branch and this the non branching side
                    
                    extend((direction ? frontSequence : backSequence), direction);
//                            std::cout<<"sequence: "<<(!side ? frontSequence : backSequence)<<std::endl;
            
                    
                }else if (edgeCounts[direction] > 1){ // if branch, we keep track of neighbours, otherwise it's a dead end and we pick another node
                    residualEdges[pair->first] = std::make_tuple(pair->second,idCounter,direction); // we preserve the edge
                }
            }
            DBGsubgraph->erase(pair->first);
//                std::cout<<*sequence->sequence<<std::endl;
        }else{
            residualEdges[pair->first] = std::make_tuple(pair->second,idCounter,0);
        }
        
        Sequence* sequence = new Sequence {std::to_string(idCounter++), "", new std::string(revCom(backSequence) + frontSequence.substr(k))}; // add sequence
        std::vector<Tag> inTags = {Tag{'f',"DP",std::to_string(pair->second.cov)},Tag{'Z',"CB",colorPalette(pair->second.color)}};
        sequence->seqPos = seqPos++; // remember the order
        GFAsubgraph.appendSegment(sequence, inTags);
    }
    jobWait(threadPool);
    while (residualEdges.size() != 0) { // construct the edges
        
        auto pair = *residualEdges.begin();
        std::string thisSegmentHeader = std::to_string(std::get<1>(pair.second));
        
        for (uint8_t i = 0; i<4; ++i) { // forward edges
            if (std::get<0>(pair.second).fw[i] != 0) {

                uint8_t nextKmer[k];
                std::string firstKmer = reverseHash(pair.first);
                firstKmer.push_back(itoc[i]);
                for (uint8_t e = 0; e<k; ++e)
                    nextKmer[e] = ctoi[(unsigned char)firstKmer[e+1]];
                
                bool isFw = false;
                uint64_t key = hash(nextKmer, &isFw);
                auto got = residualEdges.find(key);
                if (got == residualEdges.end())
                    continue;
                DBGsubgraph->erase(key);
                std::string nextSegmentHeader = std::to_string(std::get<1>(got->second));

                std::vector<Tag> inTags = {Tag{'i',"KC",std::to_string(std::get<0>(pair.second).fw[i])}};
                InEdge edge(idCounter++, edgeCounter, headersToIds[thisSegmentHeader], headersToIds[nextSegmentHeader], std::get<2>(pair.second) ? '+' : '-', std::get<2>(got->second) ? '-' : '+', std::to_string(k-1)+"M", "edge." + std::to_string(edgeCounter), inTags);
                ++edgeCounter;

                GFAsubgraph.appendEdge(edge);
            }
        }
        
        for (uint8_t i = 0; i<4; ++i) { // reverse edges
            if (std::get<0>(pair.second).bw[i] != 0) {

                uint8_t nextKmer[k];
                std::string firstKmer;
                firstKmer.push_back(itoc[i]);
                firstKmer.append(reverseHash(pair.first));
                
                for (uint8_t e = 0; e<k; ++e)
                    nextKmer[e] = ctoi[(unsigned char)firstKmer[e]];
                
                bool isFw = false;
                uint64_t key = hash(nextKmer, &isFw);
                auto got = residualEdges.find(key);
                if (got == residualEdges.end())
                    continue;
                DBGsubgraph->erase(key);
                std::string prevSegmentHeader = std::to_string(std::get<1>(got->second));

                std::vector<Tag> inTags = {Tag{'i',"KC",std::to_string(std::get<0>(pair.second).bw[i])}};
                InEdge edge(idCounter++, edgeCounter, headersToIds[prevSegmentHeader], headersToIds[thisSegmentHeader], std::get<2>(got->second) ? '+' : '-', std::get<2>(pair.second) ? '-' : '+', std::to_string(k-1)+"M", "edge." + std::to_string(edgeCounter), inTags);
                ++edgeCounter;
                GFAsubgraph.appendEdge(edge);
            }
        }
        residualEdges.erase(residualEdges.begin());
    }
    
}

void DBG::DBGgraphToGFA() {
    
    if (!userInput.noCollapse) {
        
        collapseNodes();
        
    }else{
        
        uint32_t idCounter = 0, seqPos = 0, edgeCounter = 0;
        phmap::flat_hash_map<std::string, unsigned int>& headersToIds = *GFAsubgraph.getHash1();

        phmap::parallel_flat_hash_map<uint64_t, std::string> headerLookupTable;
        
        for (auto pair : *DBGsubgraph) { // first create all nodes
            headerLookupTable[pair.first] = std::to_string(idCounter); // keep track of node ids and hashes
            std::string* inSequence = new std::string(reverseHash(pair.first));
            Sequence* sequence = new Sequence {std::to_string(idCounter++), "", inSequence};
            sequence->seqPos = seqPos++; // remember the order
            std::vector<Tag> inTags = {Tag{'f',"DP",std::to_string(pair.second.cov)},Tag{'Z',"CB",colorPalette(pair.second.color)}};
            GFAsubgraph.appendSegment(sequence, inTags);
        }
        jobWait(threadPool);
        for (auto pair : *DBGsubgraph) { // next create all edges
            
            std::string thisSegmentHeader = headerLookupTable[pair.first];
            
            for (uint8_t i = 0; i<4; ++i) { // forward edges
                if (pair.second.fw[i] != 0) {
                    
                    uint8_t nextKmer[k];
                    std::string firstKmer = reverseHash(pair.first);
                    firstKmer.push_back(itoc[i]);
                    for (uint8_t e = 0; e<k; ++e)
                        nextKmer[e] = ctoi[(unsigned char)firstKmer[e+1]];
                    
                    bool isFw = false;
                    uint64_t hashValue = hash(nextKmer, &isFw);
                    auto got = headerLookupTable.find(hashValue);
                    if (got == headerLookupTable.end())
                        continue;
                    std::string nextSegmentHeader = got->second;
                    
                    std::vector<Tag> inTags = {Tag{'i',"KC",std::to_string(pair.second.fw[i])}};
                    InEdge edge(idCounter++, edgeCounter, headersToIds[thisSegmentHeader], headersToIds[nextSegmentHeader], '+', isFw ? '+' : '-', std::to_string(k-1)+"M", "edge." + std::to_string(edgeCounter), inTags);
                    ++edgeCounter;
                    GFAsubgraph.appendEdge(edge);
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
                    
                    bool isFw = false;
                    uint64_t hashValue = hash(nextKmer, &isFw);
                    auto got = headerLookupTable.find(hashValue);
                    if (got == headerLookupTable.end())
                        continue;
                    std::string prevSegmentHeader = got->second;
                    std::vector<Tag> inTags = {Tag{'i',"KC",std::to_string(pair.second.bw[i])}};
                    InEdge edge(idCounter++, edgeCounter, headersToIds[prevSegmentHeader], headersToIds[thisSegmentHeader], isFw ? '+' : '-', '+', std::to_string(k-1)+"M", "edge." + std::to_string(edgeCounter), inTags);
                    ++edgeCounter;
                    GFAsubgraph.appendEdge(edge);
                }
            }
        }
    }
    jobWait(threadPool);
    GFAsubgraph.updateStats(); // compute summary statistics
    Report report;
    report.reportStats(GFAsubgraph, 0, 0);
}
