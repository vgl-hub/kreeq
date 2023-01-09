#include <stdlib.h>
#include <string>
#include <iostream>
#include <vector>
#include <math.h>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"

#include "kcount.h"

inline size_t Kcount::hash(uint8_t *kmer) {
    size_t result = 0;
    for(uint8_t c = 0; c<k; c++) {
        
        result += *kmer++ * pow(4,c);
    }
    return result;
}

bool Kcount::countBuff(buf64* buf, phmap::flat_hash_map<uint64_t, uint64_t>& map) {

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    for (uint64_t c = 0; c<buf->pos; c++)
        ++map[buf->seq[c]];
    
    return true;

}

bool Kcount::countUnique(phmap::flat_hash_map<uint64_t, uint64_t>& map) {
    
    uint64_t KmersUnique = 0;
    
    for (auto const &pair: map)
        ++KmersUnique;
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    totKmersUnique += KmersUnique;
    
    lck.unlock();
    
    return true;

}

void Kcount::count(std::vector<InSegment*>* segments) {
    
    uint64_t moduloMap = (uint64_t) pow(4,k) / mapCount;
    
    lg.verbose("Counting with " + std::to_string(mapCount) + " maps");
    
    for (InSegment* segment : *segments) {
        
        if (segment->getSegmentLen()<k)
            continue;
        
        uint64_t len = segment->getSegmentLen()-k+1;
        
        totKmers += len;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        uint8_t* str = new uint8_t[segment->getSegmentLen()];
        
        for (uint64_t i = 0; i < len+k-1; i++){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        for (uint64_t c = 0; c<len; ++c){
            
            uint64_t value = hash(str+c);
            
            uint16_t i = value / moduloMap;
            
            buf64* b = &buf[i];
            
            if (b->pos == b->size) {
                
                size_t newSize = b->size * 2;
                uint64_t* bufNew = new uint64_t[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            b->seq[b->pos++] = value;
            
        }
        
        delete[] str;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
    
    for(uint8_t m = 0; m<mapCount; m++)
        threadPool.queueJob([=]{ return countBuff(&buf[m], map[m]); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    for(uint8_t m = 0; m<mapCount; m++)
        threadPool.queueJob([=]{ return countUnique(map[m]); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    std::cout<<"Total kmers: "<<totKmers<<std::endl;
    std::cout<<"Unique kmers: "<<totKmersUnique<<std::endl;
    
}
