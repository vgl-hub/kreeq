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
    
    parallelMap *map;
    
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
            
            DBGkmer khmer;
            
            if (it != map->end()) {
                khmer = it->second;
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
    
    std::vector<InSegment*> inSegments = *genome->getInSegments();

    for (InSegment *inSegment : inSegments)
        jobs.push_back([this, inSegment] { return DBGtoVariants(inSegment); });

    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    deleteMapRange(mapRange);
    
}

bool DBG::searchGraph(std::array<uint16_t, 2> mapRange) { // stub
    
    parallelMap* genomeDBG = new parallelMap;
    
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
        parallelMap *map;
        
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
    parallelMap *map = maps[i];
    auto got = map->find(key);
    if (got != map->end())
        return std::make_pair(&(got->second), isFw);
    else
        return std::make_pair(nullptr,false);
    
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
                    
                    uint8_t nextKmer[k];
                    memcpy(nextKmer, origin+1, k-1);
                    nextKmer[k-1] = i;

                    if(depth == 3 && i == *(target+1)) {
                        threadLog.add("found INS\t" + newPath.sequence);
                        DBGpaths.push_back(DBGpath(INS, newPath.pos, newPath.sequence.substr(0, newPath.sequence.size()-1)));
                        break;
                    }
                    
                    if(depth == 2 && i == *target) {
                        threadLog.add("found DEL\t" + newPath.sequence);
                        newPath.type = DEL;
                        DBGpaths.push_back(newPath);
                    }
                    
                    newPath.sequence+=itoc[i];
                    std::deque<DBGpath> newDBGpaths = findPaths(nextKmer, target, depth-1, newPath, threadLog);
                    DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                    if (DBGpaths.size() > 2) // limit the number of paths to avoid extensive search
                        return DBGpaths;
                    
                    if (i == *(target+1) && newPath.sequence.size() > 1 && depth == 2) {
                        threadLog.add("found SNV\t" + newPath.sequence);
                        DBGpaths.push_back(DBGpath(SNV, newPath.pos, newPath.sequence.substr(0, newPath.sequence.size()-1)));
                    }
                }
            }
        }
    }
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

bool DBG::DBGtoVariants(InSegment *inSegment) {
    
    Log threadLog;
    threadLog.setId(inSegment->getuId());
        
    std::string sHeader = inSegment->getSeqHeader();
    parallelMap *map;
    uint64_t key, i;
    bool isFw = false;
    std::vector<std::deque<DBGpath>> variants;
        
    uint64_t len = inSegment->getSegmentLen();
    
    if (len<k)
        return true;
    
    uint64_t kcount = len-k+1;
    
    unsigned char* first = (unsigned char*)inSegment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];
    
    for (uint64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    std::vector<std::vector<uint8_t>> altPaths;
    StringGraph stringGraph(str, k);
    
    while(stringGraph.currentPos() < kcount){
        
        std::deque<DBGpath> DBGpaths;
        stringGraph.appendNext();
        altPaths = stringGraph.walkStringGraph(stringGraph.root, std::vector<uint8_t>());
//                        printAltPaths(altPaths, threadLog);
        
        bool backtrack = true;
        
        for (std::vector<uint8_t> altPath : altPaths) {
            key = hash(&altPath[0], &isFw);
            i = key % mapCount;

            map = maps[i];
            auto got = map->find(key);
                
            // check for DBG consistency
            if (got != map->end()) {
                DBGkmer &dbgkmer = got->second;
                if (stringGraph.currentPos() < kcount-1 && ((isFw && dbgkmer.fw[altPath[k]] == 0) || (!isFw && dbgkmer.bw[3-altPath[k]] == 0))) { // find alternative paths
                    std::deque<DBGpath> newDBGpaths = findPaths(&altPath[0], &altPath[k], 3, DBGpath(stringGraph.currentPos()), threadLog);
                    DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                    threadLog.add("Found " + std::to_string(DBGpaths.size()) + " alternative paths");
                    if (DBGpaths.size() > 1) { // only attempt to correct unique paths
                        DBGpaths.clear();
                        break;
                    }
                }else{backtrack = false;}
            }else{backtrack = false;}
        }
        
        uint8_t backtrackCnt = 0;
        
        if (DBGpaths.size() == 0 && backtrack) { // backtrack
            
            for (uint8_t b = 0; b < userInput.depth; ++b) {
                
                ++backtrackCnt;
                
                if (variants.size() == 0 || stringGraph.currentPos()-backtrackCnt > variants.back()[0].pos) { // prevent going past the previously discovered variant
                    
                    threadLog.add("Anomaly detected but no path is found. Backtracking at:\t" + sHeader + "\t" + std::to_string(stringGraph.currentPos()-1));
                    stringGraph.backtrack(str, k, 1);
                    altPaths = stringGraph.walkStringGraph(stringGraph.root, std::vector<uint8_t>());
                    std::vector<uint8_t> altPath = altPaths[0];
                    //                                printAltPaths(altPaths, threadLog);
                    std::deque<DBGpath> newDBGpaths = findPaths(&altPath[0], &altPath[k], 3, DBGpath(stringGraph.currentPos()), threadLog);
                    DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                    if (DBGpaths.size() > 0)
                        break;
                }
            }
        }
    
        if (DBGpaths.size() != 0) {

            // create edge at error in GFA
            threadLog.add("Candidate error at:\t" + sHeader + "\t" + std::to_string(stringGraph.currentPos()));
            std::vector<uint8_t> alts;
            
            for (DBGpath dbgpath : DBGpaths) {
                
                if (dbgpath.type == SNV) {
                    alts.push_back(stringGraph.peek());
                    alts.push_back(ctoi[(unsigned char)dbgpath.sequence[0]]);
                }else if (dbgpath.type == INS) {
                    alts.push_back(stringGraph.peek());
                    alts.push_back(4);
                }else if (dbgpath.type == DEL) {
                    alts.push_back(ctoi[(unsigned char)dbgpath.sequence[0]]);
                }
            }
            stringGraph.appendAlts(alts);
            variants.push_back(DBGpaths);
        }
        if (backtrack)
            stringGraph.advancePos(backtrackCnt);
        stringGraph.pop_front();

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
