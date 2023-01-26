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

    hashSequences(readBatch);
    
    delete readBatch;
    
    return true;
    
}

void Kpos::appendSequence(Sequence* sequence) { // method to append a new sequence from a fasta
        
    threadPool.queueJob([=]{ return inSequences.traverseInSequence(sequence); });
    
    if(verbose_flag) {std::cerr<<"\n";};
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    lck.unlock();
    
}

void Kpos::validate(UserInputKreeq userInput) {
    
    lg.verbose("Navigating with " + std::to_string(mapCount) + " maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return joinBuff(m); });
    
    jobWait(threadPool);
    
    lg.verbose("Generate summary statistics");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return stats(map[m]); });
    
    jobWait(threadPool);
    
}

bool Kpos::joinBuff(uint16_t m) {
    
//    buf64* thisBuf;
//    
//    phmap::flat_hash_map<uint64_t, std::vector<pos>>* thisMap;
//    
//    for(buf64* buf : buffers) {
//        
//        thisBuf = &buf[m];
//        
//        thisMap = &map[m];
//        
//        uint64_t len = thisBuf->pos;
//        
//        for (uint64_t c = 0; c<len; ++c)
//            ++(*thisMap)[thisBuf->seq[c]];
//        
//        delete[] thisBuf->seq;
//        
//    }
    
    return true;
    
    
}

bool Kpos::stats(phmap::flat_hash_map<uint64_t, std::vector<pos>>& map) {
    
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

void Kpos::hashSequences(Sequences* readBatch) {
    
    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    buf64* buf = new buf64[mapCount];
    
    for (Sequence* sequence : readBatch->sequences) {
        
        uint64_t len = sequence->sequence->size(), kcount = len-k+1;
        
        if (len<k)
            continue;
        
        unsigned char* first = (unsigned char*)sequence->sequence->c_str();
        
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        uint64_t value, i, newSize;
        buf64* b;
        uint64_t* bufNew;
        
        for (uint64_t c = 0; c<kcount; ++c){
            
            value = hash(str+c);
            
            i = value / moduloMap;
            
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new uint64_t[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            b->seq[b->pos++] = value;
                        
        }
        
        delete[] str;
        
        threadLog.add("Processed sequence: " + sequence->header);
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
    logs.push_back(threadLog);
    
}
