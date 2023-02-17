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

double kmerQV(uint64_t errorKmers, uint64_t totalKmers, uint8_t k){
    
    return -10*log10(1 - pow(1 - (double) errorKmers/totalKmers, (double) 1/k));
    
}

bool Kpos::traverseInReads(Sequences* readBatch) { // traverse the read

//    hashSequences(readBatch);
//
    delete readBatch;

    return true;
    
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
        uint64_t key, i, newSize;
        Buf<Gkmer>* b;
        Gkmer* bufNew;
        
        for (uint64_t c = 0; c<kcount; ++c){
            
            key = hash(str+c);
            i = key / moduloMap;
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
            gkmer->hash = key;
            gkmer->sUId = segment->getuId();
                        
        }
        
        delete[] str;
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
}

bool DBG::traverseInReads(Sequences* readBatch) { // traverse the read

    hashSequences(readBatch);

    delete readBatch;

    return true;
    
}

void DBG::hashSequences(Sequences* readBatch) {
    
    Log threadLog;
    
    threadLog.setId(readBatch->batchN);
    
    Buf<DBGkmer>* buf = new Buf<DBGkmer>[mapCount];
    
    for (Sequence* sequence : readBatch->sequences) {
        
        uint64_t len = sequence->sequence->size(), kcount = len-k+1;
        
        if (len<k)
            continue;
        
        unsigned char* first = (unsigned char*)sequence->sequence->c_str();
        
        uint8_t* str = new uint8_t[len];
        
        for (uint64_t i = 0; i<len; ++i){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        DBGkmer* dbgkmer;
        uint64_t key, i, newSize;
        Buf<DBGkmer>* b;
        DBGkmer* bufNew;
        bool isFw = false;
        bool *isFwPtr = &isFw;
        
        for (uint64_t c = 0; c<kcount; ++c){
            
            key = hash(str+c, isFwPtr);
            
            i = key / moduloMap;
            
            b = &buf[i];
            
            if (b->pos == b->size) {
                
                newSize = b->size * 2;
                bufNew = new DBGkmer[newSize];

                memcpy(bufNew, b->seq, b->size*sizeof(uint64_t));

                b->size = newSize;
                delete [] b->seq;
                b->seq = bufNew;

            }
            
            dbgkmer = &b->seq[b->pos++];
            
            if (isFwPtr)
                dbgkmer->fw[*(str+c+k)] = true;
            else
                dbgkmer->bw[*(str+c+k)] = true;
            
            dbgkmer->hash = key;
                        
        }
        
        delete[] str;
        
        threadLog.add("Processed sequence: " + sequence->header);
        
    }
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
    logs.push_back(threadLog);
    
}

void DBG::build() {
    
    lg.verbose("Navigating with " + std::to_string(mapCount) + " maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return joinBuff(m); });
    
    jobWait(threadPool);
    
    lg.verbose("Computing summary statistics");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return histogram(m); });
    
    jobWait(threadPool);
    
    uint64_t missing = pow(4,k)-totKmersDistinct;
    
    std::cout<<"DBG Summary statistics:\n"
             <<"Total: "<<totKmers<<"\n"
             <<"Unique: "<<totKmersUnique<<"\n"
             <<"Distinct: "<<totKmersDistinct<<"\n"
             <<"Missing: "<<missing<<"\n";
    
}

bool DBG::joinBuff(uint16_t m) {
    
    Buf<DBGkmer>* thisBuf;
    phmap::flat_hash_map<uint64_t, DBGkmer>* thisMap;
    
    for(Buf<DBGkmer>* buf : buffers) {
        
        thisBuf = &buf[m];
        thisMap = &map[m];
        uint64_t len = thisBuf->pos;
        
        DBGkmer* dbgkmer;
        
        for (uint64_t c = 0; c<len; ++c) {
            
            dbgkmer = &(*thisMap)[thisBuf->seq[c].hash];
            
            ++dbgkmer->cov;
            
        }
        
        delete[] thisBuf->seq;
        
    }
    
    return true;
    
    
}

bool DBG::histogram(uint16_t m) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0;
    phmap::flat_hash_map<uint64_t, uint64_t> hist;
    
    for (auto pair : map[m]) {
        
        if (pair.second.cov == 1)
            ++kmersUnique;
        
        ++kmersDistinct;
        ++hist[pair.second.cov];
        
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

void DBG::validateSequences(InSequences& inSequences) {
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    
    for (InSegment* segment : *segments) {
        
        threadPool.queueJob([=]{ return validateSegment(segment); });
        
        std::unique_lock<std::mutex> lck(mtx);
        for (auto it = logs.begin(); it != logs.end(); it++) {
         
            it->print();
            logs.erase(it--);
            
        }
        
    }
    
    jobWait(threadPool);
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        
    }
    
    std::cout<<"Presence QV (k="<<std::to_string(k)<<")\n"
             <<totErrorKmers<<"\t"
             <<totKcount<<"\t"
             <<kmerQV(totErrorKmers, totKcount, k)
             <<std::endl;
    
}

bool DBG::validateSegment(InSegment* segment) {
    
    Log threadLog;
    
    std::vector<uint64_t> errorKmers;
    int64_t len = segment->getSegmentLen(), kcount = len-k+1;
    
    if (kcount<1)
        return true;
    
    unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];
    
    for (int64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    uint64_t key, i;
    
    for (int64_t c = 0; c<kcount; ++c){
        
        key = hash(str+c);
        i = key / moduloMap;
        
        if (map[i][key].cov == 0) {
            errorKmers.push_back(c);
            threadLog.add(segment->getInSequence().substr(c, k) + " is an invalid kmer. Cov is: " + std::to_string(map[i][key].cov));
        }
    
    }
    
    threadLog.add("Processed segment: " + segment->getSeqHeader());
    threadLog.add("Found " + std::to_string(errorKmers.size()) + " error kmers out of " + std::to_string(kcount) + " kmers (presence QV: " + std::to_string(kmerQV(errorKmers.size(), kcount, k)) + ")");
    
    delete[] str;
    
    std::unique_lock<std::mutex> lck(mtx);
    
    totErrorKmers += errorKmers.size();
    totKcount += kcount;
    
    logs.push_back(threadLog);
    
    return true;
    
}
