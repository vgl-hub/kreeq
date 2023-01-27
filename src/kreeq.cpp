#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <functional>

#include <parallel_hashmap/phmap.h>
#include "parallel_hashmap/phmap_dump.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"

#include "input.h"
#include "kmer.h"
#include "kreeq.h"

bool Kpos::traverseInReads(Sequences* readBatch) { // traverse the read

//    hashSequences(readBatch);
//
//    delete readBatch;
//
//    return true;
    
}

void Kpos::index() {
    
    lg.verbose("Navigating with " + std::to_string(mapCount) + " maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return joinBuff(m); });
    
    jobWait(threadPool);
    
    lg.verbose("Compute summary statistics");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return histogram(map[m]); });
    
    jobWait(threadPool);
    
    uint64_t missing = pow(4,k)-totKmersDistinct;
    
    lg.verbose("Total: " + std::to_string(totKmers) + "\n" +
               "Unique: " + std::to_string(totKmersUnique) + "\n" +
               "Distinct: " + std::to_string(totKmersDistinct) + "\n" +
               "Missing: " + std::to_string(missing) + "\n");
    
}

bool Kpos::joinBuff(uint16_t m) {
    
    Buf<Gkmer>* thisBuf;
    phmap::flat_hash_map<uint64_t, std::vector<Gkmer*>>* thisMap;
    
    for(Buf<Gkmer>* buf : buffers) {
        
        thisBuf = &buf[m];
        thisMap = &map[m];
        uint64_t len = thisBuf->pos;
        
        for (uint64_t c = 0; c<len; ++c)
            (*thisMap)[thisBuf->seq[c].hash].push_back(&thisBuf->seq[c]);
        
        delete[] thisBuf->seq;
        
    }
    
    return true;
    
    
}

bool Kpos::histogram(phmap::flat_hash_map<uint64_t, std::vector<Gkmer*>>& map) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0;
    phmap::flat_hash_map<uint64_t, uint64_t> hist;
    
    for (auto pair : map) {
        
        if (pair.second.size() == 1)
            ++kmersUnique;
        
        ++kmersDistinct;
        ++hist[pair.second.size()];
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    totKmersUnique += kmersUnique;
    totKmersDistinct += kmersDistinct;
    
    for (auto pair : hist) {
        
        histogram1[pair.first] += pair.second;
        totKmers += pair.first * pair.second;
        
    }
    
    return true;
    
}

void Kpos::report(UserInputKreeq userInput) {
    
    const static phmap::flat_hash_map<std::string,int> string_to_case{
        {"stats",1},

    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    lg.verbose("Writing ouput: " + ext);
    
    // here we create a smart pointer to handle any kind of output stream
    std::unique_ptr<std::ostream> ostream;
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
            
        case 1: { // .stats
            
            std::ofstream ofs(userInput.outFile);
            
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            
            *ostream<<"\nTotal: "<<totKmers<<"\n";
            *ostream<<"Unique: "<<totKmersUnique<<"\n";
            *ostream<<"Distinct: "<<totKmersDistinct<<"\n";
            uint64_t missing = pow(4,k)-totKmersDistinct;
            *ostream<<"Missing: "<<missing<<"\n";
            
            ofs.close();
            
            break;
            
        }
        default: {
            
        }
            
    }
    
}

void Kpos::hashSegments() {
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    Buf<Gkmer>* buf = new Buf<Gkmer>[mapCount];
    
    for (InSegment* segment : *segments) {
        
        uint64_t len = segment->getSegmentLen(), kcount = len-k+1;
        
        if (len<k)
            continue;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i)
            str[i] = ctoi[*(first+i)];
        
        Gkmer* gkmer;
        uint64_t value, i, newSize;
        Buf<Gkmer>* b;
        Gkmer* bufNew;
        
        for (uint64_t c = 0; c<kcount; ++c){
            
            value = hash(str+c);
            i = value / moduloMap;
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new Gkmer[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            gkmer = &b->seq[b->pos++];
            gkmer->base = first+c;
            gkmer->hash = value;
            gkmer->sUId = segment->getuId();
                        
        }
        
        std::cout<<std::endl;
        
        delete[] str;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
}
