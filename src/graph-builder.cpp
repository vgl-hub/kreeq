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

bool DBG::hashSequences() {
    //   Log threadLog;
    std::string *readBatch;
    uint64_t len;
    
    while (true) {
            
        {
            
            std::unique_lock<std::mutex> lck(readMtx);
            
            if (readingDone && readBatches.size() == 0)
                return true;

            if (readBatches.size() == 0)
                continue;
            
            std::condition_variable &mutexCondition = threadPool.getMutexCondition();
            mutexCondition.wait(lck, [] {
                return !freeMemory;
            });
            
            readBatch = readBatches.front();
            readBatches.pop();
            len = readBatch->size();
            
        }
        
        if (len<k) {
            delete readBatch;
            freed += len * sizeof(char);
            continue;
        }
        
        Buf<uint8_t> *buffers = new Buf<uint8_t>[mapCount];
        const unsigned char *first = (unsigned char*) readBatch->c_str();
        allocMemory(len * sizeof(uint8_t));
        uint8_t *str = new uint8_t[len];
        uint8_t e = 0;
        uint64_t key, pos = 0, kcount = len-k+1;
        bool isFw = false;
        Buf<uint8_t>* buffer;
        
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

            key = hash(str+p, &isFw);

            buffer = &buffers[key % mapCount];
            pos = buffer->newPos(9);
            memcpy(&buffer->seq[pos-9], &key, 8);
            edgeBit edges;
            
            if (isFw){
                if (ctoi[*(first+p+k)] <= 3)
                    edges.assign(ctoi[*(first+p+k)]);
                if (p > 0 && *(str+p-1) <= 3)
                    edges.assign(4+*(str+p-1));
            }else{
                if (p > 0 && *(str+p-1) <= 3)
                    edges.assign(3-*(str+p-1));
                if (ctoi[*(first+p+k)] <= 3)
                    edges.assign(4+3-ctoi[*(first+p+k)]);
            }
            
            memcpy(&buffer->seq[pos-1], &edges, 1);
            
        }
        
        delete[] str;
        delete readBatch;
        
        //    threadLog.add("Processed sequence: " + sequence->header);
        //    std::lock_guard<std::mutex> lck(mtx);
        //    logs.push_back(threadLog);
        
        std::lock_guard<std::mutex> lck(hashMtx);
        freed += len * sizeof(char) * 2;
        buffersVec.push_back(buffers);
        
    }
    
    return true;
    
}

bool DBG::processBuffers(uint16_t m) {
    
    uint64_t pos = 0, hash;
    Buf<uint8_t> *buf;
    edgeBit edges;
    
    std::string fl = userInput.prefix + "/.buf." + std::to_string(m) + ".bin";
    std::ifstream bufFile(fl, std::ios::in | std::ios::binary);
//    map.reserve(flSize / 17); // 8 + 8 + 1
    
    while(bufFile && !(bufFile.peek() == EOF)) {
        
        {
            std::unique_lock<std::mutex> lck(mtx);
            std::condition_variable &mutexCondition = threadPool.getMutexCondition();
            mutexCondition.wait(lck, [] {
                return !freeMemory;
            });
        }
        
        parallelMap& map = *maps[m]; // the map associated to this buffer
        parallelMap32& map32 = *maps32[m];
        uint64_t map_size = mapSize(map);
        
        bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
        
        buf = new Buf<uint8_t>(pos);
        
        buf->pos = pos;
        buf->size = pos;
        
        bufFile.read(reinterpret_cast<char *>(buf->seq), sizeof(uint8_t) * buf->pos);
        
        for (uint64_t c = 0; c<pos; c+=9) {
            
            memcpy(&hash, &buf->seq[c], 8);
            memcpy(&edges, &buf->seq[c+8], 1);
            
            DBGkmer &dbgkmer = map[hash];
            bool overflow = (dbgkmer.cov == 255 ? true : false);
            
            if (dbgkmer.cov + 1 == 255)
                overflow = true;
            
            for (uint64_t w = 0; w<4; ++w) { // check weights
            
                if (dbgkmer.fw[w] + 1 == 255 || dbgkmer.bw[w] + 1 == 255) {
                    overflow = true;
                    break;
                }
            }
            
            if (!overflow) {
                
                for (uint64_t w = 0; w<4; ++w) { // update weights
                        dbgkmer.fw[w] += edges.read(w);
                        dbgkmer.bw[w] += edges.read(4+w);
                }
                
                ++dbgkmer.cov; // increase kmer coverage
            }
            
            if (overflow) {
                
                DBGkmer32 &dbgkmer32 = map32[hash];
                
                if (dbgkmer32.cov == 0) { // first time we add the kmer
                    
                    dbgkmer32 = dbgkmer;
                    dbgkmer.cov = 255; // invalidates int8 kmer
                }
                
                for (uint64_t w = 0; w<4; ++w) { // update weights
                    
                    if (LARGEST - dbgkmer32.fw[w] >= edges.read(w))
                        dbgkmer32.fw[w] += edges.read(w);
                    if (LARGEST - dbgkmer32.bw[w] >= edges.read(4+w))
                        dbgkmer32.bw[w] += edges.read(4+w);
                }
                if (dbgkmer32.cov < LARGEST)
                    ++dbgkmer32.cov; // increase kmer coverage
            }
        }
        
        delete[] buf->seq;
        freed += buf->size * sizeof(uint8_t);
        delete buf;
        alloc += mapSize(*maps[m]) - map_size;
        
        if (freeMemory || !bufFile || bufFile.peek() == EOF) {
            dumpTmpMap(userInput.prefix, m);
            reloadMap32(m);
        }
    }

    bufFile.close();
    remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str());
    
    return true;
    
}

bool DBG::reloadMap32(uint16_t m) {
    
    parallelMap& map = *maps[m]; // the map associated to this buffer
    parallelMap32& map32 = *maps32[m];
    
    for (auto pair : map32) {
        
        DBGkmer dbgkmer;
        dbgkmer.cov = 255;
        auto newPair = std::make_pair(pair.first, dbgkmer);
        map.insert(newPair);
    }
    return true;
}

bool DBG::summary(uint16_t m) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0, edgeCount = 0;
    phmap::parallel_flat_hash_map<uint64_t, uint64_t> hist;
    
    for (auto pair : *maps[m]) {
        
        if (pair.second.cov == 255) // check the large table
            continue;
        
        if (pair.second.cov == 1)
            ++kmersUnique;
        
        for (uint64_t w = 0; w<4; ++w) // update weights
            edgeCount += pair.second.fw[w] > 0 ? 1 : 0 + pair.second.bw[w] > 0 ? 1 : 0;
        
        ++kmersDistinct;
        ++hist[pair.second.cov];
    }
    
    for (auto pair : *maps32[m]) {
        
        for (uint8_t w = 0; w<4; ++w) // update weights
            edgeCount += pair.second.fw[w] > 0 ? 1 : 0 + pair.second.bw[w] > 0 ? 1 : 0;
        
        ++kmersDistinct;
        ++hist[pair.second.cov];
    }
 
    std::lock_guard<std::mutex> lck(mtx);
    totKmersUnique += kmersUnique;
    totKmersDistinct += kmersDistinct;
    totEdgeCount += edgeCount;
    
    for (auto pair : hist) {
        
        finalHistogram[pair.first] += pair.second;
        totKmers += pair.first * pair.second;
    }
    
    return true;
    
}

void DBG::DBstats() {
    
    uint64_t missing = pow(4,k)-totKmersDistinct;
    
    std::cout<<"DBG Summary statistics:\n"
             <<"Total kmers: "<<totKmers<<"\n"
             <<"Unique kmers: "<<totKmersUnique<<"\n"
             <<"Distinct kmers: "<<totKmersDistinct<<"\n"
             <<"Missing kmers: "<<missing<<"\n"
             <<"Total edges: "<<totEdgeCount<<"\n";
    
}

void DBG::kunion(){ // concurrent merging of the maps that store the same hashes
    
    parallelMap32 map32Total; // first merge high-copy kmers
    
    for (unsigned int i = 0; i < userInput.kmerDB.size(); ++i) { // for each kmerdb loads the map and merges it
        
        std::string prefix = userInput.kmerDB[i]; // loads the next map
        prefix.append("/.map.hc.bin");
        
        parallelMap32 nextMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap.phmap_load(ar_in);
        
        for (auto pair : nextMap) {
            
            DBGkmer32& dbgkmerMap32 = map32Total[pair.first];
            
            for (uint8_t w = 0; w<4; ++w) { // update weights
                
                if (LARGEST - dbgkmerMap32.fw[w] >= pair.second.fw[w])
                    dbgkmerMap32.fw[w] += pair.second.fw[w];
                else
                    dbgkmerMap32.fw[w] = LARGEST;
                if (LARGEST - dbgkmerMap32.bw[w] >= pair.second.bw[w])
                    dbgkmerMap32.bw[w] += pair.second.bw[w];
                else
                    dbgkmerMap32.bw[w] = LARGEST;
            }
            
            if (LARGEST - dbgkmerMap32.cov >= pair.second.cov)
                dbgkmerMap32.cov += pair.second.cov; // increase kmer coverage
            else
                dbgkmerMap32.cov = LARGEST;
        }
    }
    
    for (auto pair : map32Total) {
        uint64_t i = pair.first % mapCount;
        maps32[i]->insert(pair);
    }
    
    std::vector<std::function<bool()>> jobs;
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.kmerDB[0] + "/.map." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
//    for(uint32_t i : idx)
//        mergeMaps(i);
    
    dumpHighCopyKmers();
    
}

bool DBG::mergeSubMaps(parallelMap* map1, parallelMap* map2, uint8_t subMapIndex, uint16_t m) {
    
    auto& inner = map1->get_inner(subMapIndex);   // to retrieve the submap at given index
    auto& submap1 = inner.set_;        // can be a set or a map, depending on the type of map1
    auto& inner2 = map2->get_inner(subMapIndex);
    auto& submap2 = inner2.set_;
    parallelMap32& map32 = *maps32[m];
    
    for (auto pair : submap1) { // for each element in map1, find it in map2 and increase its value
        
        bool overflow = false;
        
        if (pair.second.cov == 255) // already added to int32 map
            continue;
        
        auto got = map32.find(pair.first); // check if this is already a high-copy kmer
        if (got != map32.end()) {
            overflow = true;
        }else{
            
            auto got = submap2.find(pair.first); // insert or find this kmer in the hash table
            if (got == submap2.end()) {
                submap2.insert(pair);
            }else{
                
                DBGkmer& dbgkmerMap = got->second;
                    
                if (255 - dbgkmerMap.cov <= pair.second.cov)
                    overflow = true;
                
                for (uint8_t w = 0; w<4; ++w) { // check weights
                    
                    if (255 - dbgkmerMap.fw[w] <= pair.second.fw[w] || 255 - dbgkmerMap.bw[w] <= pair.second.bw[w]) {
                        overflow = true;
                        break;
                    }
                }
                
                if (!overflow) {
                    
                    for (uint8_t w = 0; w<4; ++w) { // update weights
                        dbgkmerMap.fw[w] += pair.second.fw[w];
                        dbgkmerMap.bw[w] += pair.second.bw[w];
                    }
                    dbgkmerMap.cov += pair.second.cov; // increase kmer coverage
                }
            }
        }
        
        if (overflow) {
            
            DBGkmer32& dbgkmerMap32 = map32[pair.first];
            
            if (dbgkmerMap32.cov == 0) { // first time we add the kmer
                auto got = submap2.find(pair.first);
                DBGkmer& dbgkmerMap = got->second;
                dbgkmerMap32 = dbgkmerMap;
                dbgkmerMap.cov = 255; // invalidates int8 kmer
            }
            
            for (uint8_t w = 0; w<4; ++w) { // update weights
                
                if (LARGEST - dbgkmerMap32.fw[w] >= pair.second.fw[w])
                    dbgkmerMap32.fw[w] += pair.second.fw[w];
                else
                    dbgkmerMap32.fw[w] = LARGEST;
                if (LARGEST - dbgkmerMap32.bw[w] >= pair.second.bw[w])
                    dbgkmerMap32.bw[w] += pair.second.bw[w];
                else
                    dbgkmerMap32.bw[w] = LARGEST;
            }
            
            if (LARGEST - dbgkmerMap32.cov >= pair.second.cov)
                dbgkmerMap32.cov += pair.second.cov; // increase kmer coverage
            else
                dbgkmerMap32.cov = LARGEST;
        }
    }
    return true;
}

// subgraph functions

bool DBG::mergeSubMaps(parallelMap32* map1, parallelMap32* map2, uint8_t subMapIndex) {
    
    auto& inner = map1->get_inner(subMapIndex);   // to retrieve the submap at given index
    auto& submap1 = inner.set_;        // can be a set or a map, depending on the type of map1
    auto& inner2 = map2->get_inner(subMapIndex);
    auto& submap2 = inner2.set_;
    
    for (auto pair : submap1) { // for each element in map1, find it in map2 and increase its value
        
        auto got = submap2.find(pair.first); // insert or find this kmer in the hash table
        if (got == submap2.end()){
            submap2.insert(pair);
        }else{

            DBGkmer32& dbgkmerMap = got->second;
        
            for (uint64_t w = 0; w<4; ++w) { // update weights
                
                if (LARGEST - dbgkmerMap.fw[w] >= pair.second.fw[w])
                    dbgkmerMap.fw[w] += pair.second.fw[w];
                else
                    dbgkmerMap.fw[w] = LARGEST;
                if (LARGEST - dbgkmerMap.bw[w] >= pair.second.bw[w])
                    dbgkmerMap.bw[w] += pair.second.bw[w];
                else
                    dbgkmerMap.bw[w] = LARGEST;
                
            }
            
            if (LARGEST - dbgkmerMap.cov >= pair.second.cov)
                dbgkmerMap.cov += pair.second.cov; // increase kmer coverage
            else
                dbgkmerMap.cov = LARGEST;
            
        };
        
    }
    
    return true;
    
}


bool DBG::unionSum(parallelMap32* map1, parallelMap32* map2) {
    
    std::vector<std::function<bool()>> jobs;
    
    if (map1->subcnt() != map2->subcnt()) {
        fprintf(stderr, "Maps don't have the same numbers of submaps (%zu != %zu). Terminating.\n", map1->subcnt(), map2->subcnt());
        exit(EXIT_FAILURE);
    }
    
    for(std::size_t subMapIndex = 0; subMapIndex < map1->subcnt(); ++subMapIndex)
        jobs.push_back([this, map1, map2, subMapIndex] { return this->mergeSubMaps(map1, map2, subMapIndex); });
    
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
    return true;
    
}
