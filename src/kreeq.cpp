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
#include "node-graph.h"

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

std::vector<DBGpath> DBG::findPaths(uint8_t *origin, uint8_t *target, uint8_t depth, DBGpath currentPath) {
    
    uint8_t breadth = 0;
    std::vector<DBGpath> DBGpaths;
    
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
                        lg.verbose("found INS\t" + newPath.sequence);
                        DBGpaths.push_back(DBGpath(INS, newPath.sequence.substr(0, newPath.sequence.size()-1)));
                        break;
                    }
                    
                    if(depth == 2 && i == *target) {
                        lg.verbose("found DEL\t" + newPath.sequence);
                        newPath.type = DEL;
                        DBGpaths.push_back(newPath);
                    }
                    
                    newPath.sequence+=itoc[i];
                    std::vector<DBGpath> newDBGpaths = findPaths(nextKmer, target, depth-1, newPath);
                    DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                    if (DBGpaths.size() > 2) // limit the number of paths to avoid extensive search
                        return DBGpaths;
                    
                    if (i == *(target+1) && newPath.sequence.size() > 1 && depth == 2) {
                        lg.verbose("found SNV\t" + newPath.sequence);
                        DBGpaths.push_back(DBGpath(SNV, newPath.sequence.substr(0, newPath.sequence.size()-1)));
                    }
                }
            }
        }
    }
    return DBGpaths;
    
}

void DBG::printAltPaths(std::vector<std::vector<uint8_t>> altPaths) {
    
    uint8_t pathCounter = 1;
    
    for (std::vector<uint8_t> altPath : altPaths) {
        
        std::cout<<"P"<<+pathCounter++<<":";
        
        for (uint8_t base : altPath)
            std::cout<<itoc[base];
        std::cout<<std::endl;
    }
}

bool DBG::DBGtoGFA(std::array<uint16_t, 2> mapRange) {
    
    parallelMap* genomeDBG = new parallelMap;
    
    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    
    for (InPath& path : inPaths) {
        
        unsigned int cUId = 0, gapLen = 0;
        
        std::vector<PathComponent> pathComponents = path.getComponents();
        
        uint64_t absPos = 0;
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            
            if (component->type == SEGMENT) {
                
                uint32_t segmentCounter = 0, edgeCounter = 0;
                uint64_t cleaved = 0;
                
                auto it = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                InSegment *inSegment = *it;
                std::string sHeader = inSegment->getSeqHeader();
                parallelMap *map;
                uint64_t key, i;
                bool isFw = false;
                
                if (component->orientation == '+') {
                    
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
                    
                    bool backtrack = false;
                    
                    while(stringGraph.currentPos() < kcount){
                        
                        std::vector<DBGpath> DBGpaths;
                        stringGraph.appendNext();
                        altPaths = stringGraph.walkStringGraph(stringGraph.root, std::vector<uint8_t>());
//                        printAltPaths(altPaths);
                        
                        bool anomaly = false;
                        
                        for (std::vector<uint8_t> altPath : altPaths) {
                            
                            key = hash(&altPath[0], &isFw);
                            i = key % mapCount;
                            
                            if (i >= mapRange[0] && i < mapRange[1]) {
                                
                                map = maps[i];
                                auto got = map->find(key);
                                    
                                // check for DBG consistency
                                if (got != map->end()) {
                                    genomeDBG->insert(*got);
                                    DBGkmer &dbgkmer = got->second;
                                    if (stringGraph.currentPos() < kcount-1 && ((isFw && dbgkmer.fw[altPath[k]] == 0) || (!isFw && dbgkmer.bw[3-altPath[k]] == 0))) { // find alternative paths
                                        std::vector<DBGpath> newDBGpaths = findPaths(&altPath[0], &altPath[k], 3, DBGpath());
                                        DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                                        lg.verbose("Found " + std::to_string(DBGpaths.size()) + " alternative paths");
                                        if (DBGpaths.size() > 1) { // only attempt to correct unique paths
                                            DBGpaths.clear();
                                            anomaly = false;
                                            break;
                                        }
                                    }else{anomaly = true;}
                                }else{anomaly = true;}
                            }
                        }
                        
                        if (DBGpaths.size() == 0 && anomaly)
                            backtrack = true;
                        
                        if (backtrack) { // backtrack
                            
                            for (uint8_t b = 0; b < 5; ++b) {
                                
                                lg.verbose("Anomaly detected but no path is found. Backtracking at:\t" + sHeader + "\t" + std::to_string(stringGraph.currentPos()));
                                stringGraph.backtrack(str, k, 1);
                                altPaths = stringGraph.walkStringGraph(stringGraph.root, std::vector<uint8_t>());
                                std::vector<uint8_t> altPath = altPaths[0];
//                                printAltPaths(altPaths);
                                std::vector<DBGpath> newDBGpaths = findPaths(&altPath[0], &altPath[k], 3, DBGpath());
                                DBGpaths.insert(DBGpaths.end(), newDBGpaths.begin(), newDBGpaths.end());
                                if (DBGpaths.size() > 0)
                                    break;
                            }
                            
                            backtrack = false;
                            
                        }
                    
                        if (DBGpaths.size() != 0) {
                            
                            // create edge at error in GFA
                            lg.verbose("Candidate error at:\t" + sHeader + "\t" + std::to_string(stringGraph.currentPos()));
                            std::string newSegment1 = sHeader + "." + std::to_string(segmentCounter++);
                            std::string newSegment2 = sHeader + "." + std::to_string(segmentCounter);
                            std::string newEdge = sHeader + ".edge." + std::to_string(edgeCounter++);
                            lg.verbose("New segments after cleaving:\t" + newSegment1 + "\t" + newSegment2);
                            std::pair<InSegment*,InSegment*> segments = genome->cleaveSegment(cUId, stringGraph.currentPos()-cleaved, newSegment1, newSegment2, newEdge);
                            
                            cleaved += segments.first->getSegmentLen();
                            lg.verbose("Cleaved bases: " + std::to_string(cleaved));
                            inSegment = segments.second;
                            cUId = inSegment->getuId();
                            
                            std::pair<InSegment*,InSegment*> segments2;
                            
                            uint8_t altCounter = 0;
                            bool isCleaved = false;
                            std::vector<uint8_t> alts;
                            
                            for (DBGpath dbgpath : DBGpaths) {
                                
                                uint32_t sUId;
                                
                                if (dbgpath.type != DEL && !isCleaved) {
                                    
                                    newSegment1 = sHeader + "." + std::to_string(segmentCounter++);
                                    newSegment2 = sHeader + "." + std::to_string(segmentCounter);
                                    newEdge = sHeader + ".edge." + std::to_string(edgeCounter++);
                                    lg.verbose("New segments after cleaving:\t" + newSegment1 + "\t" + newSegment2);
                                    segments2 = genome->cleaveSegment(cUId, stringGraph.currentPos()-cleaved+1, newSegment1, newSegment2, newEdge);
                                    
                                    cleaved += segments2.first->getSegmentLen();
                                    inSegment = segments2.second;
                                    cUId = inSegment->getuId();
                                    isCleaved = true;
                                    
                                }
                                
                                if (dbgpath.type == SNV || dbgpath.type == DEL) {
                                    
                                    sUId = genome->uId.get();
                                    std::string* inSequence = new std::string(dbgpath.sequence);
                                    Sequence* sequence = new Sequence{sHeader + "." + std::to_string(segmentCounter) + ".alt" + std::to_string(++altCounter), "Candidate sequence", inSequence, NULL};
                                    lg.verbose("New segment from SNV/DEL error:\t" + sequence->header + "\t" + newSegment2);
                                    genome->traverseInSegment(sequence, std::vector<Tag>());
                                }
                                
                                if (dbgpath.type == SNV) {
                                    
                                    alts.push_back(stringGraph.peek());
                                    alts.push_back(ctoi[(unsigned char)dbgpath.sequence[0]]);
                                    
                                    InEdge edge;
                                    edge.newEdge(genome->uId.next(), segments.first->getuId(), sUId, '+', '+', "0M", sHeader + ".edge." + std::to_string(edgeCounter++));
                                    genome->appendEdge(edge);
                                    
                                    edge.newEdge(genome->uId.next(), sUId, segments2.second->getuId(), '+', '+', "0M", sHeader + ".edge." + std::to_string(edgeCounter++));
                                    genome->appendEdge(edge);
                                    
                                }else if (dbgpath.type == INS) {
                                    
                                    alts.push_back(stringGraph.peek());
                                    alts.push_back(4);
                                    
                                    InEdge edge;
                                    edge.newEdge(genome->uId.next(), segments.first->getuId(), segments2.second->getuId(), '+', '+', "0M", sHeader + ".edge." + std::to_string(edgeCounter++));
                                    genome->appendEdge(edge);
                                    
                                }else if (dbgpath.type == DEL) {
                                    
                                    alts.push_back(ctoi[(unsigned char)dbgpath.sequence[0]]);
                                    
                                    InEdge edge;
                                    edge.newEdge(genome->uId.next(), segments.first->getuId(), sUId, '+', '+', "0M", sHeader + ".edge." + std::to_string(edgeCounter++));
                                    genome->appendEdge(edge);
                                    
                                    edge.newEdge(genome->uId.next(), sUId, segments.second->getuId(), '+', '+', "0M", sHeader + ".edge." + std::to_string(edgeCounter++));
                                    genome->appendEdge(edge);
                                    
                                }
                            }
                            stringGraph.appendAlts(alts);
                        }
                        ++absPos;
                        stringGraph.pop_front();
                    }
                    delete[] str;
                }else{} // GFA not handled yet
            }else if (component->type == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                gapLen += inGap->getDist(component->start - component->end);
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
        }
    }
    delete genomeDBG;
    return true;
}
