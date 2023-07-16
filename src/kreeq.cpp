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

void DBG::initHashing(){
    
    dumpMaps = false;
    readingDone = false;
    
    for (uint8_t t = 0; t < 2; t++) {
        uint32_t jid = threadPool.queueJob([=]{ return hashSequences(); });
        dependencies.push_back(jid);
    }
    
    uint8_t threadN = threadPool.totalThreads() - 4;
    double mapsN = pow(10,log10(mapCount)/threadN), t = 0;
    
    std::array<uint16_t, 2> mapRange = {0,0};
    
    while(mapRange[1] < mapCount) {
                
        mapRange[0] = mapRange[1];
        mapRange[1] = std::ceil(pow(mapsN,t++));
        
        if (mapRange[0] >= mapRange[1])
            mapRange[1] = mapRange[0] + 1;
        
        if (mapRange[1] >= mapCount)
            mapRange[1] = mapCount;
        
        uint32_t jid = threadPool.queueJob([=]{ return processBuffers(mapRange); });
        dependencies.push_back(jid);
        
    }
    
}

bool DBG::traverseInReads(std::string* readBatch) { // specialized for string objects
    
    {
        std::lock_guard<std::mutex> lck(mtx);
        readBatches.push(readBatch);
        
        std::cout<<readBatches.size()<<std::endl;
    }
    
    mutexBatches.notify_one();
    
    return true;
    
}

bool DBG::hashSequences() {
    //   Log threadLog;
    
    std::string *readBatch;
    
    while (true) {
            
        {
            
            std::unique_lock<std::mutex> lck(mtx);
            
            if (readingDone && readBatches.size() == 0)
                return true;

            mutexBatches.wait(lck, [this] {
                return !readBatches.empty();
            });
            
            readBatch = readBatches.front();
            readBatches.pop();
            
        }
        
        uint64_t len = readBatch->size();
        
        if (len<k) {
            delete readBatch;
            return true;
        }
        
        unsigned char *first = (unsigned char*) readBatch->c_str();
        uint8_t *str = new uint8_t[len];
        uint8_t e = 0;
        uint64_t kcount = len-k+1;
        bool isFw = false;
        Buf<kmer> *buf = new Buf<kmer>(kcount);
        
        for (uint64_t p = 0; p<kcount; ++p) {
            
            for (uint8_t c = e; c<k; ++c) { // generate k bases if e=0 or the next if e=k-1
                
                str[p+c] = ctoi[*(first+p+c)]; // convert the next base
                if (str[p+c] > 3) { // if non-canonical base is found
                    p = p+c; // move position
                    e = 0; // reset base counter
                    break;
                }
                
                e = k-1;
                
            }
            
            if (e == 0) // not enough bases for a kmer
                continue;
            
            kmer &khmer = buf->seq[buf->pos++];
            
            khmer.hash = hash(str+p, &isFw);
            
            if (isFw){
                if (ctoi[*(first+p+k)] <= 3)
                    khmer.fw[ctoi[*(first+p+k)]] = 1;
                if (p > 0 && *(str+p-1) <= 3)
                    khmer.bw[*(str+p-1)] = 1;
            }else{
                if (p > 0 && *(str+p-1) <= 3)
                    khmer.fw[3-*(str+p-1)] = 1;
                if (ctoi[*(first+p+k)] <= 3)
                    khmer.bw[3-ctoi[*(first+p+k)]] = 1;
            }
            
        }
        
        delete[] str;
        delete readBatch;
        
        //    threadLog.add("Processed sequence: " + sequence->header);
        //    std::lock_guard<std::mutex> lck(mtx);
        //    logs.push_back(threadLog);
        
        std::lock_guard<std::mutex> lck(mtx);
        alloc += buf->size * sizeof(kmer);
        buffers.push_back(buf);
        
    }
    
    return true;
    
}

bool DBG::processBuffers(std::array<uint16_t, 2> mapRange) {
    
    uint16_t i;
    uint32_t b = 0;
    int64_t initial_size = 0, final_size = 0;
    Buf<kmer>* buf;
    
    while (true) {
        
        if (dumpMaps) {
            
            for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
                updateMap(userInput.prefix, m);
            
            return true;
            
        }
        
        {
            
            std::lock_guard<std::mutex> lck(mtx);
            
            alloc += final_size - initial_size;
            initial_size = 0, final_size = 0;
            
            if (readingDone && buffers.size() == 0)
                return true;
            
            if(b >= buffers.size())
                continue;
            
            buf = buffers[b];
            ++b;
            
        }
        
        for (uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
            initial_size += mapSize(*maps[m]);
        
        uint64_t len = buf->pos; // how many positions in the buffer have data
        
        for (uint64_t c = 0; c<len; ++c) {
            
            kmer &khmer = buf->seq[c];
            
            i = khmer.hash / moduloMap;
            
            if (i >= mapRange[0] && i < mapRange[1]) {
                
                phmap::flat_hash_map<uint64_t, DBGkmer>& thisMap = *maps[i]; // the map associated to this buffer
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
            
        }
        
        for (uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
            final_size += mapSize(*maps[m]);
        
    }
    
    return true;
    
}

void DBG::cleanup() {
    
    if(tmp && userInput.inDBG != userInput.prefix) {
        
        lg.verbose("Deleting tmp files");
        
        for(uint16_t m = 0; m<mapCount; ++m) // remove tmp files
            threadPool.queueJob([=]{ return remove((userInput.prefix + "./.kmap." + std::to_string(m) + ".bin").c_str()); });
        
        jobWait(threadPool);
        
        if (userInput.prefix != ".")
            rm_dir(userInput.prefix.c_str());
        
    }
    
}

void DBG::consolidate() {
    
    threadPool.status();
    
    if (!memoryOk()) {
        
        dumpMaps = true;
        readingDone = true;
        
        jobWait(threadPool, dependencies);
        
        for (Buf<kmer> *buffer : buffers) {
            
            freed += buffer->size * sizeof(kmer);
            delete[] buffer->seq;
            delete buffer;
            
        }
        
        buffers.clear();
        
        tmp = true;
        
        initHashing();
        
    }

}

bool DBG::updateMap(std::string prefix, uint16_t m) {
    
    uint64_t map_size = mapSize(*maps[m]);
    
    prefix.append("/.kmap." + std::to_string(m) + ".bin");
    
    if (fileExists(prefix)) {
    
        phmap::flat_hash_map<uint64_t, DBGkmer> *dumpMap = new phmap::flat_hash_map<uint64_t, DBGkmer>;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        dumpMap->phmap_load(ar_in);
    
        unionSum(*maps[m], *dumpMap); // merges the current map and the existing map
    
        phmap::BinaryOutputArchive ar_out(prefix.c_str()); // dumps the data
        dumpMap->phmap_dump(ar_out);
    
        delete dumpMap;
    
    }else{
    
        phmap::BinaryOutputArchive ar_out(prefix.c_str()); // dumps the data
        maps[m]->phmap_dump(ar_out);
    
    }
    
    freed += map_size;
    delete maps[m];
    maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
    
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

void DBG::summary() {
    
    {
        std::lock_guard<std::mutex> lck(mtx);
        readingDone = true;
    }
    
    jobWait(threadPool, dependencies);
    
    for (Buf<kmer> *buffer : buffers) {
        
        freed += buffer->size * sizeof(kmer);
        delete[] buffer->seq;
        delete buffer;
        
    }
    
    if (tmp) {
        
        for(uint16_t m = 0; m<mapCount; ++m) {
            uint32_t jid = threadPool.queueJob([=]{ return updateMap(userInput.prefix, m); });
            dependencies.push_back(jid);
        }
        
        lg.verbose("Updating maps");
        
        jobWait(threadPool, dependencies);
        
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
            
            for(uint16_t m = mapRange[0]; m<=mapRange[1]; ++m)
                threadPool.queueJob([=]{ return loadMap(userInput.prefix, m); });
            
            jobWait(threadPool);
            
        }
        
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

bool DBG::validateSegment(InSegment* segment, std::array<uint16_t, 2> mapRange) {
    
    std::vector<uint64_t> missingKmers;
    std::vector<uint64_t> edgeMissingKmers;
    
    uint64_t len = segment->getSegmentLen();
    
    if (len<k)
        return true;
    
    uint64_t kcount = len-k+1, kmers = 0;
    
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
        
//        std::cout<<"\n"<<itoc[*(str+c)]<<"\t"<<c<<"\t"<<isFw<<std::endl;
        
        if (i >= mapRange[0] && i <= mapRange[1]) {
            
            map = maps[i];
            
            auto it = map->find(key);
            const DBGkmer *dbgkmer = (it == map->end() ? NULL : &it->second);
            
            if (it == map->end()) // merqury QV
                missingKmers.push_back(c);
            else if (it->second.cov < userInput.covCutOff) // merqury QV with cutoff
                missingKmers.push_back(c);
            else if (dbgkmer != NULL) { // kreeq QV
                
                if (isFw){
                    
                    if ((c<kcount-1 && dbgkmer->fw[*(str+c+k)] == 0) && (c>0 && dbgkmer->bw[*(str+c-1)] == 0)){
                        edgeMissingKmers.push_back(c);
//                        std::cout<<"edge error1"<<std::endl;
                    }
                }else{
                    
                    if ((c>0 && dbgkmer->fw[3-*(str+c-1)] == 0) && (c<kcount-1 && dbgkmer->bw[3-*(str+c+k)] == 0)){
                        edgeMissingKmers.push_back(c);
//                        std::cout<<"edge error2"<<std::endl;
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
    totEdgeMissingKmers += edgeMissingKmers.size();
    
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

void DBG::load() {
    
    tmp = true;
    userInput.prefix = userInput.inDBG;
    
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
    
    lg.verbose("Writing ouput: " + userInput.outFile);
    
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
