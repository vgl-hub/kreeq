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

uint64_t mapSize(phmap::flat_hash_map<uint64_t, DBGkmer>& m) {
    
    return (m.size() * (sizeof(DBGkmer) + sizeof(void*)) + // data list
     m.bucket_count() * (sizeof(void*) + sizeof(uint64_t))) // bucket index
    * 1.5; // estimated allocation overheads
    
}

double errorRate(uint64_t missingKmers, uint64_t totalKmers, uint8_t k){ // estimate QV from missing kmers
    
    return 1 - pow(1 - (double) missingKmers/totalKmers, (double) 1/k);
    
}

bool DBG::memoryOk() {
    
    return get_mem_inuse(3) < (userInput.maxMem == 0 ? get_mem_total(3) * 0.5 : userInput.maxMem);
    
}

bool DBG::memoryOk(int64_t delta) {
    
    return get_mem_inuse(3) + convert_memory(delta, 3) < (userInput.maxMem == 0 ? get_mem_total(3) * 0.5 : userInput.maxMem);
    
}

bool DBG::traverseInReads(std::string* readBatch) { // specialized for string objects
    
    alloc += readBatch->size() * sizeof(char);
    
    uint32_t jid = threadPool.queueJob([=]{ return hashSequences(readBatch); });
    dependencies.push_back(jid);
    
    return true;
    
}

bool DBG::hashSequences(std::string* readBatch) {
    
    Log threadLog;
    
    Buf<kmer>* buf = new Buf<kmer>[mapCount];
        
    uint64_t len = readBatch->size();
    
    if (len<k)
        return true;
    
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
            
            kmer* khmer;
            uint64_t key, i, newSize;
            Buf<kmer>* b;
            kmer* bufNew;
            bool isFw = false;
            
            for (uint64_t c = 0; c<kcount; ++c){
                
                key = hash(str+c+p-e, &isFw);
                i = key / moduloMap;
                b = &buf[i];
                
                if (b->pos == b->size) {
                    
                    newSize = b->size*2;
                    bufNew = new kmer[newSize];
                    
                    memcpy(bufNew, b->seq, b->size*sizeof(kmer));
                    
                    b->size = newSize;
                    delete[] b->seq;
                    b->seq = bufNew;
                    
                }
                
                khmer = &b->seq[b->pos++];
                
                if (isFw){
                    khmer->fw[*(str+c+k+p-e)] = 1;
                    if (c > 0)
                        khmer->bw[*(str+c-1+p-e)] = 1;
                }else{
                    if (c > 0)
                        khmer->fw[3-*(str+c-1+p-e)] = 1;
                    khmer->bw[3-*(str+c+k+p-e)] = 1;
                }
                
                khmer->hash = key;
                
            }
            
            e = 0;
            
        }else{
            
            ++e;
            
        }
        
    }

    delete[] str;
    delete readBatch;
    
    // track memory usage
    uint64_t thisAlloc = 0;
    for(uint64_t i = 0 ; i < mapCount ; ++i)
        thisAlloc += buf[i].size * sizeof(kmer);
        
    // threadLog.add("Processed sequence: " + sequence->header);
    
    freed += len * sizeof(char);
    alloc += thisAlloc;
    
    std::lock_guard<std::mutex> lck(mtx);
    buffers.push_back(buf);
    logs.push_back(threadLog);
    
    return true;
    
}

void DBG::finalize() {
    
    lg.verbose("Finalizing DBG");
    
    if (tmp) {
        
        lg.verbose("Using tmp");
        
        updateDBG();
        
        lg.verbose("DBG updated");
        
    }else{
        
        for(uint16_t m = 0; m<mapCount; ++m)
            threadPool.queueJob([=]{ return countBuffs(m); });
        
        jobWait(threadPool);
        
    }
    
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

void DBG::cleanup() {
    
    if(tmp) {
        
        lg.verbose("Deleting tmp files");
        
        for(uint16_t m = 0; m<mapCount; ++m) // remove tmp files
            threadPool.queueJob([=]{ return remove(("./.kmap." + std::to_string(m) + ".bin").c_str()); });
        
        jobWait(threadPool);
        
    }
    
}

void DBG::consolidate() { // to reduce memory footprint we consolidate the buffers as we go
    
//    for (unsigned int i = 0; i<buffers.size(); ++i) { // for each buffer, check if we can delete them
//        
//        unsigned int counter = 0;
//        
//        for(uint16_t m = 0; m<mapCount; ++m) { // for each map
//            
//            Buf<kmer> *thisBuf = &buffers[i][m];
//            
//            if(thisBuf->seq == NULL){
//                
//                ++counter; // keeps track of the buffers that were processed so far
//                
//                if (counter == mapCount) {
//                    
//                    delete[] buffers[i];
//                    buffers.erase(buffers.begin() + i);
//                    
//                }
//                
//            }
//            
//        }
//        
//    }
    
    for(uint16_t m = 0; m<mapCount; ++m) { // for each map, consolidate
     
        if (mapsInUse[m] == false) {
            
            uint32_t jid = threadPool.queueJob([=]{ return countBuffs(m); });
            dependencies.push_back(jid);
            
        }
        
    }
    
    threadPool.status();
    
    if (!memoryOk()) {
        
        updateDBG();
        tmp = true;
        
    }

}

void DBG::updateDBG() {
    
    lg.verbose("\nCompleting residual jobs");
    
    jobWait(threadPool, dependencies);
    
    for(uint16_t m = 0; m<mapCount; ++m) {
        uint32_t jid = threadPool.queueJob([=]{ return countBuffs(m); });
        dependencies.push_back(jid);
    }
    
    lg.verbose("Counting all residual buffers");
    
    jobWait(threadPool, dependencies);
    
    lg.verbose("Removing residual heap memory allocations");
    
    for(Buf<kmer>* buf : buffers)
        delete[] buf;
    
    buffers.clear();
    
    for(uint16_t m = 0; m<mapCount; ++m) {
        uint32_t jid = threadPool.queueJob([=]{ return updateMap(userInput.prefix, m); });
        dependencies.push_back(jid);
    }
    
    lg.verbose("Updating maps");
    
    jobWait(threadPool, dependencies);
    
}

bool DBG::updateMap(std::string prefix, uint16_t m) {
    
    uint64_t map_size1 = 0, map_size2 = 0;
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::flat_hash_map<uint64_t, DBGkmer> dumpMap;
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    dumpMap.phmap_load(ar_in);
    
    map_size1 = mapSize(dumpMap);
    map_size2 = mapSize(*maps[m]);
    alloc += map_size1;
    
    unionSum(*maps[m], dumpMap); // merges the current map and the existing map
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str()); // dumps the data
    dumpMap.phmap_dump(ar_out);
    
    delete maps[m];
    maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
    
    freed += map_size1 + map_size2;
    
    return true;
    
}

bool DBG::unionSum(phmap::flat_hash_map<uint64_t, DBGkmer>& map1, phmap::flat_hash_map<uint64_t, DBGkmer>& map2) {
    
    for (auto pair : map1) { // for each element in map1, find it in map2 and increase its value
        
        DBGkmer &dbgkmerMap = map2[pair.first]; // insert or find this kmer in the hash table
        
        for (uint64_t w = 0; w<4; ++w) { // update weights
            
            if (255 - dbgkmerMap.fw[w] >= pair.second.fw[w])
                dbgkmerMap.fw[w] += pair.second.fw[w];
            else
                dbgkmerMap.fw[w] = 255;
            if (255 - dbgkmerMap.bw[w] >= pair.second.bw[w])
                dbgkmerMap.bw[w] += pair.second.bw[w];
            else
                dbgkmerMap.bw[w] = 255;
            
        }
        
        if (255 - dbgkmerMap.cov >= pair.second.cov)
            dbgkmerMap.cov += pair.second.cov; // increase kmer coverage
        else
            dbgkmerMap.cov = 255;
        
    }
    
    return true;
    
}


bool DBG::countBuffs(uint16_t m) { // counts all residual buffers for a certain map as we finalize the kmerdb

    if (mapsInUse[m] == true) {return true;}
    
    uint64_t releasedMem = 0, initial_size = 0, final_size = 0;
    
    initial_size = mapSize(*maps[m]);
    
    {
        std::lock_guard<std::mutex> lck(mtx);
        mapsInUse[m] = true;
    }

    for(Buf<kmer>* buf : buffers) {
                    
        releasedMem = buf[m].size * sizeof(kmer);
        
        countBuff(&buf[m], m);
        
    }
    
    final_size = mapSize(*maps[m]);
    
    alloc += final_size - initial_size;
    freed += releasedMem;
    std::lock_guard<std::mutex> lck(mtx); // release the map
    mapsInUse[m] = false;
    
    return true;

}

bool DBG::countBuff(Buf<kmer>* buf, uint16_t m) { // counts a single buffer
    
    Buf<kmer> &thisBuf = *buf;
    
    if (thisBuf.seq != NULL) { // sanity check that this buffer was not already processed
        
        phmap::flat_hash_map<uint64_t, DBGkmer>& thisMap = *maps[m]; // the map associated to this buffer
        
        uint64_t len = thisBuf.pos; // how many positions in the buffer have data
        
        for (uint64_t c = 0; c<len; ++c) {
            
            kmer &khmer = thisBuf.seq[c];
            DBGkmer &dbgkmer = thisMap[khmer.hash]; // insert or find this kmer in the hash table
            
            for (uint64_t w = 0; w<4; ++w) { // update weights

                if (255 - dbgkmer.fw[w] >= khmer.fw[w])
                    dbgkmer.fw[w] += khmer.fw[w];
                if (255 - dbgkmer.bw[w] >= khmer.bw[w])
                    dbgkmer.bw[w] += khmer.bw[w];
            }
            if (dbgkmer.cov < 255)
                ++dbgkmer.cov; // increase kmer coverage
            
        }
        
        delete[] thisBuf.seq; // delete the buffer
        thisBuf.seq = NULL; // set its sequence to the null pointer in case its checked again
        
    }
    
    return true;

}

bool DBG::histogram(uint16_t m) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0, map_size = 0;
    phmap::flat_hash_map<uint64_t, uint64_t> hist;
    
    if (tmp)
        loadMap(userInput.prefix, m);
    
    for (auto pair : *maps[m]) {
        
        if (pair.second.cov == 1)
            ++kmersUnique;
        
        ++kmersDistinct;
        ++hist[pair.second.cov];
        
    }
    
    if (tmp) {
        map_size = mapSize(*maps[m]);
        delete maps[m];
        maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
        freed += map_size;
    }
    
    std::lock_guard<std::mutex> lck(mtx);
    totKmersUnique += kmersUnique;
    totKmersDistinct += kmersDistinct;
    
    for (auto pair : hist) {
        
        finalHistogram[pair.first] += pair.second;
        totKmers += pair.first * pair.second;
        
    }
    
    return true;
    
}

void DBG::validateSequences(InSequences &inSequences) {
    
    lg.verbose("Validating sequence");
    
    std::vector<InSegment*> *segments = inSequences.getInSegments();
    
    std::array<uint16_t, 2> mapRange = {0,0};
    
    while(mapRange[1] < mapCount-1) {
        
        uint64_t max = 0;
        
        for(uint16_t m = mapRange[0]; m<mapCount; ++m) {
            
            if(tmp)
                max += fileSize(userInput.prefix + "/.kmap." + std::to_string(m) + ".bin");
            
            if(!memoryOk(max))
                break;
            
            mapRange[1] = m;
            
        }
    
        if(tmp) {
            
            for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
                threadPool.queueJob([=]{ return loadMap(userInput.prefix, m); });
            
            jobWait(threadPool);
            
        }
        
        std::cout<<mapRange[0]<<" "<<mapRange[1]<<std::endl;
        
        for (InSegment* segment : *segments)
            threadPool.queueJob([=]{ return validateSegment(segment, mapRange); });
        
        jobWait(threadPool);
        
        if (tmp) {
            
            for(uint16_t m = mapRange[0]; m<=mapRange[1]; ++m) {
                
                uint64_t map_size = mapSize(*maps[m]);
                delete maps[m];
                maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
                freed += map_size;
                
            }
            
            mapRange[0] = mapRange[1] + 1;
            
        }
        
    }
    
    for (auto it = logs.begin(); it != logs.end(); it++) {
     
        it->print();
        logs.erase(it--);
        
    }
    
    double merquryError = errorRate(totMissingKmers, totKcount, k), merquryQV = -10*log10(merquryError);
    
    std::cout<<"Presence QV (k="<<std::to_string(k)<<")\n"
             <<totMissingKmers<<"\t"
             <<totKcount<<"\t"
             <<merquryQV<<"\t"
             <<merquryError<<std::endl;
    
}

bool DBG::validateSegment(InSegment* segment, std::array<uint16_t, 2> mapRange) {
    
    std::vector<uint64_t> missingKmers;
    std::vector<uint64_t> edgeMissingKmers;
    uint64_t len = segment->getSegmentLen(), kcount = len-k+1, kmers = 0;
    
    if (kcount<1)
        return true;
    
    unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
    uint8_t* str = new uint8_t[len];    
    
    for (uint64_t i = 0; i<len; ++i)
        str[i] = ctoi[*(first+i)];
    
    uint64_t key, i;
    
    phmap::flat_hash_map<uint64_t, DBGkmer> *map;
    
    // kreeq QV
    bool isFw = false;
    
//    std::cout<<segment->getSeqHeader()<<std::endl;
    
    for (uint64_t c = 0; c<kcount; ++c){
        
        key = hash(str+c, &isFw);
        
        i = key / moduloMap;
        
        if (i >= mapRange[0] && i <= mapRange[1]) {
            
            map = maps[i];
            
            auto it = map->find(key);
            const DBGkmer *dbgkmer = (it == map->end() ? NULL : &it->second);
            
            //std::cout<<"\n"<<itoc[*(str+c)]<<"\t"<<c<<"\t"<<isFw<<std::endl;
            
            if (it == map->end()) // merqury QV
                missingKmers.push_back(c);
            else if (it->second.cov < userInput.covCutOff) // merqury QV with cutoff
                missingKmers.push_back(c);
            else if (dbgkmer != NULL) { // kreeq QV
                
                if (isFw){
                    
                    if ((c<kcount-1 && dbgkmer->fw[*(str+c+k)]) == 0 || (c>0 && dbgkmer->bw[*(str+c-1)] == 0)){
                        edgeMissingKmers.push_back(c);
                        //std::cout<<"edge error1"<<std::endl;
                    }
                }else{
                    
                    if ((c>0 && dbgkmer->fw[3-*(str+c-1)] == 0) || (c<kcount-1 && dbgkmer->bw[3-*(str+c+k)] == 0)){
                        edgeMissingKmers.push_back(c);
                        //std::cout<<"edge error2"<<std::endl;
                    }
                }
            }
            ++kmers;
        }
    }
    
    //double merquryError = errorRate(totMissingKmers, kcount, k), merquryQV = -10*log10(merquryError);
    //double kreeqError = errorRate(totMissingKmers + edgeMissingKmers.size(), kcount, k), kreeqQV = -10*log10(kreeqError);
    
    delete[] str;

    totMissingKmers += missingKmers.size();
    totKcount += kmers;
    
    return true;
    
}

bool DBG::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    uint64_t map_size = mapSize(*maps[m]);
    
    delete maps[m];
    maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
    
    freed += map_size;
    
    return true;
    
}

void DBG::load() { // concurrent loading of existing hashmaps
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return loadMap(userInput.iDBGFileArg, m); });
    
    jobWait(threadPool);
    
}

bool DBG::loadMap(std::string prefix, uint16_t m) { // loads a specific map
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    maps[m]->phmap_load(ar_in);
    
    alloc += mapSize(*maps[m]);
    
    return true;

}

void DBG::report() { // generates the output from the program
    
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
