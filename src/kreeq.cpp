#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <functional>
#include <list>
#include <iomanip>

#include "parallel_hashmap/phmap.h"
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

double kmerQV(uint64_t missingKmers, uint64_t totalKmers, uint8_t k){ // estimate QV from missing kmers
    
    return -10*log10(1 - pow(1 - (double) missingKmers/totalKmers, (double) 1/k));
    
}

bool DBG::traverseInReads(std::string* readBatch) { // specialized for string objects

    hashSequences(readBatch);
    
    freed = readBatch->size() * sizeof(char);
    delete readBatch;
    
    return true;
    
}

void DBG::hashSequences(std::string* readBatch) {
    
    Log threadLog;
    
    Buf<DBGkmer>* buf = new Buf<DBGkmer>[mapCount];
        
    uint64_t len = readBatch->size();
    
    if (len<k)
        return;
    
    unsigned char* first = (unsigned char*) readBatch->c_str();
    
    uint8_t* str = new uint8_t[len];
    uint64_t e = 0;
    
    for (uint64_t p = 0; p<len; ++p) {
        
        str[p] = ctoi[*(first+p)];
        
        if (str[p] > 3 || p+1 == len){
            
            if (p+1 == len && str[p] < 4) { // end of sequence, adjust indexes
                ++e;
                ++p;
            }
            
            if (e < k) { // beginning/end of a sequence or kmer too short, nothing to be done
                e = 0;
                continue;
            }
            
            uint64_t kcount = e-k+1;
            
            DBGkmer* dbgkmer;
            uint64_t key, i, newSize;
            Buf<DBGkmer>* b;
            DBGkmer* bufNew;
            bool isFw = false;
            bool *isFwPtr = &isFw;
            
            for (uint64_t c = 0; c<kcount; ++c){
                
                key = hash(str+c+p-e);
                i = key / moduloMap;
                b = &buf[i];
                
                if (b->pos == b->size) {
                    
                    newSize = b->size * 2;
                    bufNew = new DBGkmer[newSize];
                    
                    memcpy(bufNew, b->seq, b->size*sizeof(DBGkmer));
                    
                    b->size = newSize;
                    delete[] b->seq;
                    b->seq = bufNew;
                    
                }
                
                dbgkmer = &b->seq[b->pos++];
                
                if (isFwPtr)
                    dbgkmer->fw[*(str+c+k)] = 1;
                else
                    dbgkmer->bw[*(str+c+k)] = 1;
                
                dbgkmer->hash = key;
                
            }
            
            e = 0;
            
        }else{
            
            ++e;
            
        }
        
    }
    
    delete[] str;
        
//        threadLog.add("Processed sequence: " + sequence->header);
    
    std::unique_lock<std::mutex> lck(mtx);
    
    buffers.push_back(buf);
    
    logs.push_back(threadLog);
    
}

void DBG::finalize() {
    
    lg.verbose("Navigating with " + std::to_string(mapCount) + " maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return countBuffs(m); });
    
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

void DBG::consolidate() { // to reduce memory footprint we consolidate the buffers as we go
    
    for (unsigned int i = 0; i<buffers.size(); ++i) { // for each buffer
        
        unsigned int counter = 0;
        
        for(uint16_t m = 0; m<mapCount; ++m) { // for each map
            
            Buf<DBGkmer> *thisBuf = &buffers[i][m];
            
            if (thisBuf->seq != NULL && mapsInUse[m] == false) { // if the buffer was not counted and the associated map is not in use we process it
                
                mapsInUse[m] = true;
                threadPool.queueJob([=]{ return countBuff(thisBuf, m); });
                
            }
            
            if(thisBuf->seq == NULL){
                
                ++counter; // keeps track of the buffers that were processed so far
                
                if (counter == mapCount) {
                    lg.verbose("Jobs waiting/running: " + std::to_string(threadPool.queueSize()) + "/" + std::to_string(threadPool.running()) + " memory used/total: " + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
                    buffers.erase(buffers.begin() + i);
                }
                
            }

        }
        
    }

}

bool DBG::countBuffs(uint16_t m) { // counts all residual buffers for a certain map as we finalize the kmerdb
    
    for(Buf<DBGkmer>* buf : buffers)
        countBuff(&buf[m], m);
    
    return true;

}

bool DBG::countBuff(Buf<DBGkmer>* buf, uint16_t m) { // counts a single buffer
    
    Buf<DBGkmer> &thisBuf = *buf;
    
    if (thisBuf.seq != NULL) { // sanity check that this buffer was not already processed
        
        phmap::flat_hash_map<uint64_t, DBGkmer>& thisMap = map[m]; // the map associated to this buffer
        
        uint64_t len = thisBuf.pos; // how many positions in the buffer have data
        
        for (uint64_t c = 0; c<len; ++c) {
            
            DBGkmer &dbgkmerBuf = thisBuf.seq[c];
            DBGkmer &dbgkmerMap = thisMap[dbgkmerBuf.hash]; // insert or find this kmer in the hash table
            
            for (uint64_t w = 0; w<4; ++w) { // update weights
                
                if (255 - dbgkmerMap.fw[w] >= dbgkmerBuf.fw[w])
                    dbgkmerMap.fw[w] += dbgkmerBuf.fw[w];
                if (255 - dbgkmerMap.bw[w] >= dbgkmerBuf.bw[w])
                    dbgkmerMap.bw[w] += dbgkmerBuf.bw[w];
            }
            
            ++dbgkmerMap.cov; // increase kmer coverage
            
        }
        
        delete[] thisBuf.seq; // delete the buffer
        thisBuf.seq = NULL; // set its sequence to the null pointer in case its checked again
        
    }
    
    std::unique_lock<std::mutex> lck(mtx); // release the map
    mapsInUse[m] = false;
    
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
        
        finalHistogram[pair.first] += pair.second;
        totKmers += pair.first * pair.second;
        
    }
    
    return true;
    
}

void DBG::validateSequences(InSequences& inSequences) {
    
    lg.verbose("Validating sequence");
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    
    for (InSegment* segment : *segments) {
        
        //threadPool.queueJob([=]{ return
        
        validateSegment(segment);
            
        //});
        
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
             <<totMissingKmers<<"\t"
             <<totKcount<<"\t"
             <<kmerQV(totMissingKmers, totKcount, k)
             <<std::endl;
    
}

bool DBG::validateSegment(InSegment* segment) {
    
    Log threadLog;
    
    std::vector<uint64_t> missingKmers;
    int64_t len = segment->getSegmentLen(), kcount = len-k+1;
    
    if (kcount<1)
        return true;
    
    unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];
    
    for (int64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    uint64_t key, i;
    
    // kreeq QV
    double presenceProbs[] = {0.000001, 0.01, 0.1};
    double kmerProb = presenceProbs[0];
    std::list<double> kmerProbs, weightsProbs;
    std::vector<double> totProbs;
    
    std::cout<<segment->getSeqHeader()<<std::endl;
    
    for (int64_t c = 0; c<kcount; ++c){
        
        key = hash(str+c);
        i = key / moduloMap;
        
        auto it = map[i].find(key);
        
        if (it == map[i].end()) // merqury QV
            missingKmers.push_back(c);
        
        // kreeq QV
        double totProb = 0;
        DBGkmer dbgkmer = it->second;
        kmerProb = dbgkmer.cov < 3 ? presenceProbs[dbgkmer.cov] : 1;
        std::cout<<"kcov: "<<std::to_string(dbgkmer.cov)<<" "<<std::setprecision(15)<<"prob: "<<kmerProb<<"\t";
        
        kmerProbs.push_back(kmerProb);
        if (c >= k)
            kmerProbs.pop_front();
        
        for (double n : kmerProbs)
            totProb += n;
            
        for (double n : weightsProbs)
            totProb *= n;
        
        totProbs.push_back(totProb/kmerProbs.size());
        
        std::cout<<"Final probability: "<<std::setprecision(15)<<std::to_string(totProb/kmerProbs.size())<<std::endl;
    
    }
    
    threadLog.add("Processed segment: " + segment->getSeqHeader());
    threadLog.add("Found " + std::to_string(missingKmers.size()) + " missing kmers out of " + std::to_string(kcount) + " kmers (presence QV: " + std::to_string(kmerQV(missingKmers.size(), kcount, k)) + ")");
    
    delete[] str;
    
    std::unique_lock<std::mutex> lck(mtx);
    
    totMissingKmers += missingKmers.size();
    totKcount += kcount;
    
    logs.push_back(threadLog);
    
    return true;
    
}

bool DBG::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    map[m].phmap_dump(ar_out);
    
    freeContainer(map[m]);
    
    return true;
    
}

void DBG::load(UserInputKreeq& userInput) { // concurrent loading of existing hashmaps
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return loadMap(userInput.iDBGFileArg, m); });
    
    jobWait(threadPool);
    
}

bool DBG::loadMap(std::string prefix, uint16_t m) { // loads a specific maps
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    map[m].phmap_load(ar_in);
    
    histogram(m);
    
    return true;

}

void DBG::report(UserInputKreeq& userInput) { // generates the output from the program
    
    const static phmap::flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"kreeq",1}
    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    lg.verbose("Writing ouput: " + ext);
    
    std::unique_ptr<std::ostream> ostream; // smart pointer to handle any kind of output stream
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
        
        default:
        case 1: { // .kreeq
            
            make_dir(userInput.outFile.c_str());
            
            std::ofstream ofs(userInput.outFile + "/.index");
            
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            
            *ostream<<+k<<"\n"<<mapCount<<std::endl;
            
            ofs.close();
            
            for(uint16_t m = 0; m<mapCount; ++m)
                threadPool.queueJob([=]{ return dumpMap(userInput.outFile, m); }); // writes map to file concurrently
            
            jobWait(threadPool);
            
            break;
            
        }
            
    }
    
}