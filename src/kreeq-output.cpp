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

void DBG::report() { // generates the output from the program
    
    const static phmap::parallel_flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"kreeq",1},
        {"bedtable",2},
        {"csvtable",2},
        {"kwig",3},
        {"bkwig",4},
        {"gfa",5},
        {"errorbed",6}
    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    if (userInput.outFile.find(".") != std::string::npos || userInput.outFile == "" || ext == "kreeq" || userInput.stats_flag)
        DBGstats();
    
    lg.verbose("Writing ouput: " + userInput.outFile);
    
    std::unique_ptr<std::ostream> ostream; // smart pointer to handle any kind of output stream
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
        default:
        case 2:
        case 3:
        case 4: { // .bed .csv .kwig .bkwig
            validateSequences(); // validate the input sequence
            break;
        }
        case 5: {}
            
    }
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
        
        default:
        case 1: { // .kreeq
            std::ofstream ofs(userInput.outFile + "/.index"); // adding index
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            *ostream<<+k<<"\n"<<mapCount<<std::endl;
            ofs.close();
            break;
        }
        case 2: { // .bed .csv
            printTable(ext);
            break;
        }
        case 3: { // .kwig
            printTableCompressed();
            break;
        }
        case 4: { // .bkwig
            printTableCompressedBinary();
            break;
        }
        case 5: { // .gfa
            printGFA();
            break;
        }
    }
}

void DBG::printTable(std::string ext) {
    
    char colSep = ',', entrySep = ',';
    
    if (ext == "bed")
        colSep = '\t', entrySep = ':';
    else if (ext == "csv")
        colSep = ',', entrySep = ' ';
    
    std::ofstream ofs(userInput.outFile);
    
    genome->sortPathsByOriginal();
    
    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    std::vector<DBGbase*> *dbgbases = genome->getInSegmentsDBG();
    
    for (InPath& path : inPaths) {
        
        unsigned int cUId = 0, gapLen = 0, sIdx = 0;
        
        std::vector<PathComponent> pathComponents = path.getComponents();
        
        uint64_t absPos = 0;
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            
            if (component->type == SEGMENT) {
                
                auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                if (inSegment != inSegments->end()) {sIdx = std::distance(inSegments->begin(), inSegment);} // gives us the segment index
                
                DBGbase *dbgbase = (*dbgbases)[sIdx];
                
                std::vector<uint8_t> kmerCov(k-1,0);
                std::vector<uint8_t> edgeCovFw(k-1,0);
                std::vector<uint8_t> edgeCovBw(k-1,0);
                
                if (component->orientation == '+') {
                    
                    for (uint64_t i = 0; i < (*inSegment)->getSegmentLen(); ++i) {
                        
                        ofs<<path.getHeader()
                           <<colSep<<absPos<<colSep;
                        
                        kmerCov.push_back(dbgbase[i].cov);
                        
                        for(uint8_t c = 0; c<k; ++c){
                            ofs<<std::to_string(kmerCov[c]);
                            if (c < k - 1)
                                ofs<<entrySep;
                        }
                        
                        ofs<<colSep;
                        
                        edgeCovFw.push_back(dbgbase[i].isFw ? dbgbase[i].fw : dbgbase[i].bw);
                        
                        for(uint8_t c = 0; c<k; ++c){
                            ofs<<std::to_string(edgeCovFw[c]);
                            if (c < k - 1)
                                ofs<<entrySep;
                        }
                        
                        ofs<<colSep;
                        
                        edgeCovBw.push_back(dbgbase[i].isFw ? dbgbase[i].bw : dbgbase[i].fw);
                        
                        for(uint8_t c = 0; c<k; ++c){
                            ofs<<std::to_string(edgeCovBw[c]);
                            if (c < k - 1)
                                ofs<<entrySep;
                        }
                        
                        ofs<<"\n";
                        
                        kmerCov.erase(kmerCov.begin());
                        edgeCovFw.erase(edgeCovFw.begin());
                        edgeCovBw.erase(edgeCovBw.begin());
                        
                        ++absPos;
                        
                    }
                    
                }else{
                    
                    // GFA not handled yet
                    
                }
                
            }else if (component->type == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                gapLen += inGap->getDist(component->start - component->end);
                
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
        }
    }
    ofs.close();
}

void DBG::printTableCompressed() {
    
    std::ofstream ofs(userInput.outFile);
    ofs<<std::to_string(k)<<"\n";
    
    genome->sortPathsByOriginal();
    
    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    std::vector<DBGbase*> *dbgbases = genome->getInSegmentsDBG();
    
    for (InPath& path : inPaths) {
        
        unsigned int cUId = 0, gapLen = 0, sIdx = 0;
        
        std::vector<PathComponent> pathComponents = path.getComponents();
        
        uint64_t absPos = 0;
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            
            if (component->type == SEGMENT) {
                
                auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                if (inSegment != inSegments->end()) {sIdx = std::distance(inSegments->begin(), inSegment);} // gives us the segment index
                
                DBGbase *dbgbase = (*dbgbases)[sIdx];
                
                ofs<<"fixedStep chrom="<<path.getHeader()<<" start="<<absPos<<" step=1"<<"\n";
                
                if (component->orientation == '+') {
                    
                    for (uint64_t i = 0; i < (*inSegment)->getSegmentLen(); ++i) {
                        
                        ofs<<std::to_string(dbgbase[i].cov)<<","<<std::to_string(dbgbase[i].isFw ? dbgbase[i].fw : dbgbase[i].bw)<<","<<std::to_string(dbgbase[i].isFw ? dbgbase[i].bw : dbgbase[i].fw)<<"\n";
                        
                        ++absPos;
                        
                    }
                    
                }else{
                    
                    // GFA not handled yet
                    
                }
                
            }else if (component->type == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                gapLen += inGap->getDist(component->start - component->end);
                
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
        }
    }
    ofs.close();
}

void DBG::printTableCompressedBinary() {
    
    std::ofstream ofs(userInput.outFile, std::fstream::trunc | std::ios::out | std::ios::binary);
    ofs.write(reinterpret_cast<const char *>(&k), sizeof(uint8_t));
    genome->sortPathsByOriginal();
    
    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    std::vector<DBGbase*> *dbgbases = genome->getInSegmentsDBG();
    
    for (InPath& path : inPaths) {
        
        unsigned int cUId = 0, gapLen = 0, sIdx = 0;
        
        std::vector<PathComponent> pathComponents = path.getComponents();
        
        uint64_t absPos = 0;
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            
            if (component->type == SEGMENT) {
                
                auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                if (inSegment != inSegments->end()) {sIdx = std::distance(inSegments->begin(), inSegment);} // gives us the segment index
                
                DBGbase *dbgbase = (*dbgbases)[sIdx];
                
                std::string str("fixedStep chrom="+path.getHeader()+" start="+std::to_string(absPos)+" step=1");
                uint16_t size=str.size();
                ofs.write(reinterpret_cast<const char *>(&size), sizeof(uint16_t));
                ofs<<str;
                uint64_t len = (*inSegment)->getSegmentLen();
                ofs.write(reinterpret_cast<const char *>(&len), sizeof(uint64_t));
                
                if (component->orientation == '+') {
                    
                    for (uint64_t i = 0; i < len; ++i) {
                        
                        ofs.write(reinterpret_cast<const char *>(&dbgbase[i].cov), sizeof(uint8_t));
                        ofs.write(reinterpret_cast<const char *>(dbgbase[i].isFw ? &dbgbase[i].fw : &dbgbase[i].bw), sizeof(uint8_t));
                        ofs.write(reinterpret_cast<const char *>(dbgbase[i].isFw ? &dbgbase[i].bw : &dbgbase[i].fw), sizeof(uint8_t));
                        
                        ++absPos;
                        
                    }
                    
                }else{
                    
                    // GFA not handled yet
                    
                }
                
            }else if (component->type == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                gapLen += inGap->getDist(component->start - component->end);
                
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
        }
    }
    ofs.close();
}

void DBG::printGFA() {
    
    lg.verbose("Generating GFA");

//    std::vector<std::function<bool()>> jobs;
    std::array<uint16_t, 2> mapRange = {0,0};
    
    genome->sortPathsByOriginal();
    
    while (mapRange[1] < mapCount) {
        
        mapRange = computeMapRange(mapRange);
        
        loadMapRange(mapRange);
        
//        for (InPath& path : inPaths)
//            jobs.push_back([this, path, mapRange] { return DBGtoGFA(path, mapRange); });
        DBGtoGFA(mapRange);
        
//        threadPool.queueJobs(jobs);
//        jobWait(threadPool);
//        jobs.clear();
        
        deleteMapRange(mapRange);
        
    }

    Report report;
    report.outFile(*genome, userInput.outFile, userInput, 0);
    
}
