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

inline size_t Kcount::hash(unsigned short int *kmer) {
    size_t result = 0;
    for(unsigned short int c = 0; c<k; c++) {
        
        result += *kmer++ * pow(4,c);
    }
    return result;
}

bool Kcount::countBuff(std::vector<unsigned long long int>& buff, phmap::flat_hash_map<unsigned long long int, unsigned long long int>& map) {

//    only if sorted table is needed:
//    std::sort(buff.begin(), buff.end());
    
    for (unsigned long long int& value : buff)
        ++map[value];
    
    return true;

}

bool Kcount::countUnique(phmap::flat_hash_map<unsigned long long int, unsigned long long int>& map) {
    
    unsigned long long int KmersUnique = 0;
    
    for (auto const &pair: map)
        ++KmersUnique;
    
    std::unique_lock<std::mutex> lck (mtx, std::defer_lock);
    
    lck.lock();
    
    totKmersUnique += KmersUnique;
    
    lck.unlock();
    
    return true;

}

void Kcount::count(std::vector<InSegment*>* segments) {
    
    unsigned long long int moduloMap = (unsigned long long int) pow(4,k) / mapCount;
    
    lg.verbose("Counting with " + std::to_string(mapCount) + " maps");
    
    for (InSegment* segment : *segments) {
        
        if (segment->getSegmentLen()<k)
            continue;
        
        unsigned long long int len = segment->getSegmentLen()-k+1;
        
        totKmers += len;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        unsigned short int* str = new unsigned short int[segment->getSegmentLen()];
        
        for (unsigned long long int i = 0; i < len+k-1; i++){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        for (unsigned long long int c = 0; c<len; ++c){
            
            unsigned long long int value = hash(str+c);
            
            unsigned long long int i = value / moduloMap;
            
            buff[i].push_back(value);
            
        }
        
        delete[] str;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
    
    for(unsigned short int m = 0; m<mapCount; m++)
        threadPool.queueJob([=]{ return countBuff(buff[m], map[m]); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    for(unsigned short int m = 0; m<mapCount; m++)
        threadPool.queueJob([=]{ return countUnique(map[m]); });
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    std::cout<<"Total kmers: "<<totKmers<<std::endl;
    std::cout<<"Unique kmers: "<<totKmersUnique<<std::endl;
    
}
