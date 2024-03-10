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

#include "input.h"
#include "kmer.h"
#include "kreeq.h"

void DBG::report() { // generates the output from the program
    
    const static phmap::parallel_flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"kreeq",1},
        {"bed",2},
        {"csvtable",2},
        {"kwig",3},
        {"bkwig",4},
        {"gfa",5},
        {"gfa2",5},
        {"gfa.gz",5},
        {"gfa2.gz",5},
        {"vcf",6}
    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    if (userInput.outFile.find(".") != std::string::npos || userInput.outFile == "" || ext == "kreeq" || userInput.stats_flag)
        stats();
    
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
        case 5:
        case 6: { // .gfa*, .vcf
            correctSequences();
            break;
        }
            
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
        case 6: { // .vcf
            printVCF();
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
            
            if (component->componentType == SEGMENT) {
                
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
                
            }else if (component->componentType == GAP){
                
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
            
            if (component->componentType == SEGMENT) {
                
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
                
            }else if (component->componentType == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                gapLen += inGap->getDist(component->start - component->end);
                
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
        }
    }
    ofs.close();
}

void DBG::writeIndex(std::ofstream &ofs) { // writes: indexSize, nPaths, and for each path writes size and path header and nComponents, and for each segment 1) type 2) absPos 3) segment length

    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    
    uint32_t nPaths = inPaths.size();
    ofs.write(reinterpret_cast<const char *>(&nPaths), sizeof(uint32_t));
    
    for (InPath& path : inPaths) {
        
        uint32_t cUId = 0;
        std::vector<PathComponent> pathComponents = path.getComponents();
        uint64_t absPos = 0;
        
        uint16_t headerSize = path.getHeader().size();
        ofs.write(reinterpret_cast<const char *>(&headerSize), sizeof(uint16_t));
        ofs<<path.getHeader();
        
        uint32_t nComponents = pathComponents.size();
        ofs.write(reinterpret_cast<const char *>(&nComponents), sizeof(uint32_t));
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            ofs.write(reinterpret_cast<const char *>(&component->componentType), sizeof(ComponentType));
            
            if (component->componentType == SEGMENT) {
                
                auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                ofs.write(reinterpret_cast<const char *>(&absPos), sizeof(uint64_t));
                uint64_t segmentlen = (*inSegment)->getSegmentLen();
                ofs.write(reinterpret_cast<const char *>(&segmentlen), sizeof(uint64_t));
                uint8_t step = 1;
                ofs.write(reinterpret_cast<const char *>(&step), sizeof(uint8_t));
                absPos += segmentlen;
                
            }else if (component->componentType == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                ofs.write(reinterpret_cast<const char *>(&absPos), sizeof(uint64_t));
                uint64_t gapLen = inGap->getDist(component->start - component->end);
                ofs.write(reinterpret_cast<const char *>(&gapLen), sizeof(uint64_t));
                absPos += gapLen;
                
            }else{}
        }
    }
}

void DBG::printTableCompressedBinary() {
    
    std::ofstream ofs(userInput.outFile, std::fstream::trunc | std::ios::out | std::ios::binary);
    ofs.write(reinterpret_cast<const char *>(&k), sizeof(uint8_t));
    
    genome->sortPathsByOriginal();
    writeIndex(ofs);
    
    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    std::vector<DBGbase*> *dbgbases = genome->getInSegmentsDBG();
    
    for (InPath& path : inPaths) {
        
        unsigned int cUId = 0, gapLen = 0, sIdx = 0;
        std::vector<PathComponent> pathComponents = path.getComponents();
        uint64_t absPos = 0;
        
        uint16_t pHeaderLen = path.getHeader().size();
        ofs.write(reinterpret_cast<const char *>(&pHeaderLen), sizeof(uint16_t));
        ofs<<path.getHeader();
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            
            if (component->componentType == SEGMENT) {
                
                auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                if (inSegment != inSegments->end()) {sIdx = std::distance(inSegments->begin(), inSegment);} // gives us the segment index
                
                DBGbase *dbgbase = (*dbgbases)[sIdx];
                uint64_t len = (*inSegment)->getSegmentLen();
                
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
            }else if (component->componentType == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                gapLen += inGap->getDist(component->start - component->end);
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
        }
    }
    ofs.close();
}

void DBG::printGFA() {
    
    genome->sortSegmentsByOriginal();
    genome->sortEdgesByOriginal();
    genome->sortPathsByOriginal();

    Report report;
    report.outFile(*genome, userInput.outFile, userInput, 0);
    
}

void DBG::printVCF() {
    
    genome->sortPathsByOriginal();

    Report report;
    report.outFile(*genome, userInput.outFile, userInput, 0);
    
}
