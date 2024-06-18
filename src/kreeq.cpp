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

void DBG::correctSequences() {
    
    if (userInput.inSequence.empty())
        return;
    
    lg.verbose("Generating candidates for correction");
    
    std::array<uint16_t, 2> mapRange = {0,0};
    
    // this section will need to be implemented to manage memory
    
//    if (computeMapRange(mapRange)[1] < mapCount) {
//
//        for (uint8_t i = 0; i < userInput.depth; ++i) {
//
//            mapRange = {0,0};
//
//            while (mapRange[1] < mapCount) {
//
//                mapRange = computeMapRange(mapRange);
//                loadMapRange(mapRange);
//                //            searchGraph(mapRange);
//                deleteMapRange(mapRange);
//
//            }
//        }
//    }
    
    std::vector<std::function<bool()>> jobs;
    mapRange = computeMapRange(mapRange);
    loadMapRange(mapRange);
    
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
    
    std::vector<InSegment*> inSegments = *genome->getInSegments();
    for (InSegment *inSegment : inSegments)
        jobs.push_back([this, inSegment] { return DBGtoVariants(inSegment); });
    
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    deleteMapRange(mapRange);
    
}

bool DBG::searchGraph(std::array<uint16_t, 2> mapRange) { // stub
    
    ParallelMap* genomeDBG = new ParallelMap;
    
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    
    for (InSegment *inSegment : *inSegments) {
    
        uint64_t len = inSegment->getSegmentLen();
        
        if (len<k)
            return true;
        
        uint64_t kcount = len-k+1;
        
        unsigned char* first = (unsigned char*)inSegment->getInSequencePtr()->c_str();
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i)
            str[i] = ctoi[*(first+i)];
        
        uint64_t key, i;
        bool isFw = false;
        ParallelMap *map;
        
        for(uint64_t c = 0; c<kcount; ++c){
            
            key = hash(str+c, &isFw);
            i = key % mapCount;
            
            if (i >= mapRange[0] && i < mapRange[1]) {
                
                map = maps[i];
                auto got = map->find(key);
                
                // check for DBG consistency
                if (got == map->end()) {
                    genomeDBG->insert(*got);
//                    DBGkmer &dbgkmer = got->second;
                }
            }
        }
    }
    delete genomeDBG;
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

int8_t DBG::scorePath(std::string path) { // not tested
    
    int8_t score = 0;
    uint8_t len = path.size();
    unsigned char* first = (unsigned char*)path.c_str();
    
    for (uint8_t i = 0; i < len-k; ++i) {
        
//        std::cout<<itoc[first[i]];
        
        auto it = findDBGkmer(first+i);
        DBGkmer *dbgOrigin = it.first;
        
        if (dbgOrigin == nullptr) {
            --score;
        }else{
            ++score;
        }
    }
    
//    std::cout<<std::endl;
    
    return score;
    
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

std::deque<DBGpath> DBG::findPaths(uint8_t *origin, uint8_t *target, uint8_t depth, DBGpath currentPath, Log &threadLog) {
    
    uint8_t breadth = 0;
    std::deque<DBGpath> DBGpaths;
    
    for (uint8_t a = 0; a <= breadth; ++a) { // to expand the search before and after the current pos
        
        origin -= a;
        
        if (depth != 0) {
            
            auto it = findDBGkmer(origin);
            DBGkmer *dbgOrigin = it.first;
            bool isFw = it.second;
            
            if (dbgOrigin == nullptr)
                return DBGpaths;
            
            for (uint8_t i = 0; i < 4 ; ++i){
                
                DBGpath newPath = currentPath;
                
                uint8_t e = 0;
                if(isFw)
                    e = i;
                else
                    e = 3-i;
                
                if ((isFw && dbgOrigin->fw[e] != 0) || (!isFw && dbgOrigin->bw[e] != 0)) {
                    
//                    threadLog.add(std::to_string(i) +":"+ std::to_string(depth) + ":" + std::to_string(*(target+1)));
                    
                    if(depth == 3 && i == *(target+1)) {
                        newPath.score += checkNext(origin, target+1);
                            threadLog.add("found INS (score: " + std::to_string(newPath.score) + ")\t" + newPath.sequence);
                            DBGpaths.push_back(DBGpath(INS, newPath.pos, newPath.sequence.substr(0, newPath.sequence.size()-1), newPath.score));
                        newPath = currentPath;
                    }
                    
                    uint8_t nextKmer[k];
                    memcpy(nextKmer, origin+1, k-1);
                    nextKmer[k-1] = i;
                    
                    if(depth == 2 && i == *target) {
                        newPath.score += checkNext(nextKmer, target+1);
                            threadLog.add("found DEL (score: " + std::to_string(newPath.score) + ")\t" + newPath.sequence);
                            newPath.type = DEL;
                            DBGpaths.push_back(newPath);
                        newPath = currentPath;
                    }
                    
                    if (depth == 2 && i == *(target+1)) {
                        newPath.score += checkNext(nextKmer, target+2);
                            threadLog.add("found SNV (score: " + std::to_string(newPath.score) + ")\t" + newPath.sequence);
                            DBGpaths.push_back(DBGpath(SNV, newPath.pos, newPath.sequence.substr(0, newPath.sequence.size()), newPath.score));
                        newPath = currentPath;
                    }
                    
                    newPath.sequence+=itoc[i];
                    
//                    if (DBGpaths.size() > 2) // limit the number of paths to avoid extensive search
//                        return DBGpaths;
                    
                    std::deque<DBGpath> newDBGpaths = findPaths(nextKmer, target, depth-1, newPath, threadLog);
                    DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                }
            }
        }
    }
//    if (DBGpaths.size() > 2) // limit the number of paths to avoid extensive search
//        DBGpaths.clear();
    return DBGpaths;
    
}

void DBG::printAltPaths(std::vector<std::vector<uint8_t>> altPaths, Log &threadLog) {
    
    uint8_t pathCounter = 1;
    
    for (std::vector<uint8_t> altPath : altPaths) {
        
        std::string logMsg = "P" + std::to_string(pathCounter++) + ":";
        
        for (uint8_t base : altPath)
            logMsg += itoc[base];
        threadLog.add(logMsg += "\n");
    }
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

bool DBG::detectAnomalies(InSegment *inSegment, std::vector<uint64_t> &anomalies) {
    
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

bool DBG::DBGtoVariants(InSegment *inSegment) {
    
    Log threadLog;
    threadLog.setId(inSegment->getuId());
        
    std::string sHeader = inSegment->getSeqHeader();
    ParallelMap *map;
    uint64_t key, i;
    bool isFw = false;
    std::vector<std::deque<DBGpath>> variants;
    std::vector<uint64_t> anomalies;
    userInput.inBedInclude == "" ? detectAnomalies(inSegment, anomalies) : loadSegmentCoordinates(inSegment, anomalies);
    uint64_t len = inSegment->getSegmentLen();
    
    if (len<k || anomalies.size() == 0)
        return true;

    uint64_t kcount = len-k+1;
    
    unsigned char* first = (unsigned char*)inSegment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];
    
    for (uint64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    std::vector<std::vector<uint8_t>> altPaths;
    
    StringGraph *stringGraph = new StringGraph(str, k, anomalies[0]-k);
    
    for(std::vector<uint64_t>::iterator anomaly = anomalies.begin(); anomaly < anomalies.end(); ++anomaly) {
        
        std::deque<DBGpath> DBGpaths;
        altPaths = stringGraph->walkStringGraph(stringGraph->root, std::vector<uint8_t>());
//        printAltPaths(altPaths, threadLog);
        
        if (altPaths.size() < 2) {
            stringGraph->deleteStringGraph(stringGraph->root);
//                delete stringGraph;
            stringGraph = new StringGraph(str, k, *anomaly-k);
            altPaths = stringGraph->walkStringGraph(stringGraph->root, std::vector<uint8_t>());
        }
        
        bool backtrack = true;
        
        for (std::vector<uint8_t> altPath : altPaths) {
            key = hash(&altPath[0], &isFw);
            i = key % mapCount;
            
            map = maps[i];
            auto got = map->find(key);
            
            // check for DBG consistency
            if (got != map->end()) {
                
                DBGkmer &dbgkmer = got->second;
                if (stringGraph->currentPos() < kcount-1 && ((isFw && dbgkmer.fw[altPath[k]] == 0) || (!isFw && dbgkmer.bw[3-altPath[k]] == 0))) { // find alternative paths
                    
                    double score = - checkNext(&altPath[0], &altPath[k]);
                    
                    threadLog.add(std::to_string(altPath[0]) + "," + std::to_string(altPath[k]) + "," + std::to_string(score));
                    
                    std::deque<DBGpath> newDBGpaths = findPaths(&altPath[0], &altPath[k], userInput.depth, DBGpath(stringGraph->currentPos(), score), threadLog);
                    
                    DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                    threadLog.add("Found " + std::to_string(DBGpaths.size()) + " alternative paths");
                    //                    if (DBGpaths.size() > 2) { // only attempt to correct unique paths
                    //                        DBGpaths.clear();
                    //                        break;
                    //                    }
                    
                }else{backtrack = false;}
            }else{backtrack = false;}
        }
        
        if (DBGpaths.size() > 1) {
            std::sort(DBGpaths.begin(), DBGpaths.end(), [](const DBGpath& v1, const DBGpath& v2) {return v1.score > v2.score;});
            DBGpaths = std::deque<DBGpath>(DBGpaths.begin(), DBGpaths.begin()+1); // get at most two branches in the search tree
        }
        
        uint8_t backtrackCnt = 0;
        
        if (DBGpaths.size() == 0 && backtrack) { // backtrack
            
            threadLog.add("Anomaly detected but no path is found. Backtracking");
            
            for (uint8_t b = 0; b < userInput.backtrackingSpan; ++b) {
                
                ++backtrackCnt;
                
                if (variants.size() == 0 || stringGraph->currentPos()-backtrackCnt > variants.back()[0].pos) { // prevent going past the previously discovered variant
                    
                    threadLog.add("Testing position:\t" + sHeader + "\t" + std::to_string(stringGraph->currentPos()));
                    stringGraph->backtrack(str, k, 1);
                    altPaths = stringGraph->walkStringGraph(stringGraph->root, std::vector<uint8_t>());
                    std::vector<uint8_t> altPath = altPaths[0];
                    //                                printAltPaths(altPaths, threadLog);
                    double score = - checkNext(&altPath[0], &altPath[k]);
                    
                    std::deque<DBGpath> newDBGpaths = findPaths(&altPath[0], &altPath[k], userInput.depth, DBGpath(stringGraph->currentPos(), score), threadLog);
                    
                    DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                }
            }
            
            if (DBGpaths.size() > 1) {
                std::sort(DBGpaths.begin(), DBGpaths.end(), [](const DBGpath& v1, const DBGpath& v2) {return v1.score > v2.score;});
                DBGpaths = std::deque<DBGpath>(DBGpaths.begin(), DBGpaths.begin()+1); // get at most two branches in the search tree
            }
        }
        
        if (DBGpaths.size() != 0) {
            
            threadLog.add("Candidate error at:\t" + sHeader + "\t" + std::to_string(stringGraph->currentPos()+1));
            std::vector<uint8_t> alts;
            
            for (DBGpath dbgpath : DBGpaths) {
                
                if (dbgpath.type == SNV) {
                    alts.push_back(stringGraph->peek());
                    alts.push_back(ctoi[(unsigned char)dbgpath.sequence[0]]);
                }else if (dbgpath.type == INS) {
                    alts.push_back(stringGraph->peek());
                    alts.push_back(4);
                }else if (dbgpath.type == DEL) {
                    alts.push_back(ctoi[(unsigned char)dbgpath.sequence[0]]);
                }
            }
            stringGraph->appendAlts(alts);
            variants.push_back(DBGpaths);
        }
        if (backtrack)
            stringGraph->advancePos(backtrackCnt);
        if (std::next(anomaly) != anomalies.end()) {
            if (*std::next(anomaly) < (*anomaly)+k) {
                stringGraph->appendNext(*std::next(anomaly) - *anomaly);
                for (uint8_t e = 0; e < *std::next(anomaly) - *anomaly; ++e)
                    stringGraph->pop_front();
            }
        }else{
            delete stringGraph;
        }
    }
    delete[] str;
    
    inSegment->addVariants(variants);
    
    std::string ext = getFileExt("." + userInput.outFile);
    if (ext == "gfa" || ext == "gfa2" || ext == "gfa.gz" || ext == "gfa2.gz")
        variantsToGFA(inSegment, threadLog);
    
    std::unique_lock<std::mutex> lck (mtx);
    logs.push_back(threadLog);
    
    return true;
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
    mergeSubgraphs();
    DFS();
    summary(*DBGsubgraph);
    DBGgraphToGFA();
    
}

void DBG::summary(ParallelMap32& DBGsubgraph) {
    
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
    ParallelMap32 *segmentSubmap = new ParallelMap32;
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
                        segmentSubmap->insert(*got);
                    }else{
                        map32 = maps32[i];
                        auto got = map32->find(key);
                        segmentSubmap->insert(*got);
                    }
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

void DBG::DFS() {
    
    ParallelMap32 candidates, newCandidates;
    ParallelMap32* subgraph = DBGsubgraph;
    
    std::array<uint16_t, 2> mapRange = {0,0};
    for (uint8_t i = 0; i < userInput.kmerDepth; ++i) {

        mapRange = {0,0};

        while (mapRange[1] < mapCount) {

            mapRange = computeMapRange(mapRange);
            loadMapRange(mapRange);
            newCandidates = DFSpass(subgraph, mapRange);
            deleteMapRange(mapRange);
            candidates.insert(newCandidates.begin(), newCandidates.end());
            subgraph = &newCandidates;
        }
    }
    DBGsubgraph->insert(candidates.begin(), candidates.end());
}

ParallelMap32 DBG::DFSpass(ParallelMap32* subgraph, std::array<uint16_t, 2> mapRange) {
    
    ParallelMap32 newCandidates;
    
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
                
                auto got = DBGsubgraph->find(key);
                
                if (got == DBGsubgraph->end()) {
                    
                    uint64_t m = key % mapCount;
                    
                    if (m >= mapRange[0] && m < mapRange[1]) {
                        
                        map = maps[m];
                        auto got = map->find(key);
                        
                        if (got != map->end()) {
                            
                            if (got->second.cov != 255) {
                                newCandidates.insert(*got);
                            }else{
                                map32 = maps32[m];
                                auto got = map32->find(key);
                                newCandidates.insert(*got);
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
                
                auto got = DBGsubgraph->find(key);
                
                if (got == DBGsubgraph->end()) {
                    
                    uint64_t m = key % mapCount;
                    
                    if (m >= mapRange[0] && m < mapRange[1]) {
                        
                        map = maps[m];
                        auto got = map->find(key);
                        
                        if (got != map->end()) {
                            if (got->second.cov != 255) {
                                newCandidates.insert(*got);
                            }else{
                                map32 = maps32[m];
                                auto got = map32->find(key);
                                newCandidates.insert(*got);
                            }
                        }
                    }
                }
            }
        }
    }
    return newCandidates;
}

void DBG::mergeSubgraphs() {
    
    for (ParallelMap32 *map1 : DBGTmpSubgraphs) {
        unionSum(map1, DBGsubgraph);
        delete map1;
    }
    
}

void DBG::DBGgraphToGFA() {
    
    phmap::parallel_flat_hash_map<uint64_t, uint32_t> idLookupTable;
    uint32_t idCounter = 0, seqPos = 0, edgeCounter = 0;
    
    for (auto pair : *DBGsubgraph) { // first create all nodes
        idLookupTable[pair.first] = idCounter; // keep track of node ids and hashes
        std::string* inSequence = new std::string(reverseHash(pair.first));
        Sequence* sequence = new Sequence {std::to_string(idCounter++), "", inSequence};
        sequence->seqPos = seqPos; // remember the order
        std::vector<Tag> inTags = {Tag{'i',"RC",std::to_string(pair.second.cov*k)}};
        GFAsubgraph.appendSegment(sequence, inTags);
        seqPos++;
    }
    
    for (auto pair : *DBGsubgraph) { // next create all edges
    
        uint32_t thisSegmentId = idLookupTable[pair.first];
        
        for (uint8_t i = 0; i<4; ++i) { // forward edges
            if (pair.second.fw[i] != 0) {
                
                uint8_t nextKmer[k];
                std::string firstKmer = reverseHash(pair.first);
                firstKmer.push_back(itoc[i]);
                for (uint8_t e = 0; e<k; ++e)
                    nextKmer[e] = ctoi[(unsigned char)firstKmer[e+1]];

                bool isFw = false;
                uint64_t hashValue = hash(nextKmer, &isFw);
                auto got = idLookupTable.find(hashValue);
                if (got == idLookupTable.end())
                    continue;
                uint32_t nextSegmentId = got->second;
                std::vector<Tag> inTags = {Tag{'i',"KC",std::to_string(pair.second.fw[i])}};
                InEdge edge(idCounter++, edgeCounter, thisSegmentId, nextSegmentId, '+', isFw ? '+' : '-', "1N"+std::to_string(k-1)+"M", "edge." + std::to_string(edgeCounter), inTags);
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
                auto got = idLookupTable.find(hashValue);
                if (got == idLookupTable.end())
                    continue;
                uint32_t prevSegmentId = got->second;
                std::vector<Tag> inTags = {Tag{'i',"KC",std::to_string(pair.second.bw[i])}};
                InEdge edge(idCounter++, edgeCounter, prevSegmentId, thisSegmentId, isFw ? '+' : '-', '+', "1N"+std::to_string(k-1)+"M", "edge." + std::to_string(edgeCounter), inTags);
                ++edgeCounter;
                GFAsubgraph.appendEdge(edge);
            }
        }
    }
    jobWait(threadPool);
    GFAsubgraph.updateStats(); // compute summary statistics
    Report report;
    report.reportStats(GFAsubgraph, 0, 0);
}
