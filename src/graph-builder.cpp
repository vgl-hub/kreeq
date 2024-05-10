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

uint64_t mapSize(parallelMap& m) {
    
   return m.capacity() * (sizeof(parallelMap::value_type) + 1) + sizeof(parallelMap);
    
}

void DBG::status() {
    
    std::chrono::duration<double> elapsed = std::chrono::high_resolution_clock::now() - past;
    
    if (elapsed.count() > 0.1) {
        lg.verbose("Read batches: " + std::to_string(readBatches.size()) + ". Hash buffers: " + std::to_string(buffersVec.size()) + ". Memory in use/allocated/total: " + std::to_string(get_mem_inuse(3)) + "/" + std::to_string(get_mem_usage(3)) + "/" + std::to_string(get_mem_total(3)) + " " + memUnit[3], true);
    
        past = std::chrono::high_resolution_clock::now();
    }
    
}

void DBG::joinThreads() {
    
    uint8_t threadsDone = 0;
    bool done = false;
    
    while (!done) {
        for (uint8_t i = 0; i < futures.size(); ++i) {
            if (futures[i].wait_for(std::chrono::milliseconds(0)) == std::future_status::ready) {
                if (threads[i].joinable()) {
                    threads[i].join();
                    ++threadsDone;
                }
            }else{
                status();
            }
        }
        if (threadsDone == futures.size())
            done = true;
    }
    
    futures.clear();
    threads.clear();
    
}

bool DBG::memoryOk() {
    
    return get_mem_inuse(3) < maxMem;
    
}

bool DBG::memoryOk(int64_t delta) {
    
    return get_mem_inuse(3) + convert_memory(delta, 3) < maxMem;
    
}

bool DBG::traverseInReads(std::string* readBatch) { // specialized for string objects
    
    while(freeMemory) {status();}
    
    {
        std::lock_guard<std::mutex> lck(readMtx);
        readBatches.push(readBatch);
    }
    
    return true;
    
}

void DBG::consolidate() {
    
    status();

}

void DBG::initHashing(){
    
    std::packaged_task<bool()> task([this] { return dumpBuffers(); });
    futures.push_back(task.get_future());
    threads.push_back(std::thread(std::move(task)));
    
    int16_t threadN = threadPool.totalThreads() - 1; // substract the writing thread
    
    if (threadN == 0)
        threadN = 1;
    
    for (uint8_t t = 0; t < threadN; t++) {
        std::packaged_task<bool()> task([this] { return hashSequences(); });
        futures.push_back(task.get_future());
        threads.push_back(std::thread(std::move(task)));
    }
    
}

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
            pos = buffer->newPos(9);;
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

bool DBG::dumpBuffers() {
    
    bool hashing = true;
    std::vector<Buf<uint8_t>*> buffersVecCpy;
    std::ofstream bufFile[mapCount];
    
    for (uint16_t b = 0; b<mapCount; ++b) // we open all files at once
        bufFile[b] = std::ofstream(userInput.prefix + "/.buf." + std::to_string(b) + ".bin", std::fstream::app | std::ios::out | std::ios::binary);
    
    while (hashing) {
        
        if (readingDone) { // if we have finished reading the input
            
            uint8_t hashingDone = 0;
        
            for (uint8_t i = 1; i < futures.size(); ++i) { // we check how many hashing threads are still running
                if (futures[i].wait_for(std::chrono::milliseconds(0)) == std::future_status::ready)
                    ++hashingDone;
            }
            
            if (hashingDone == futures.size() - 1) // if all hashing threads are done we exit the loop after one more iteration
                hashing = false;
                    
        }
        
        {
            std::unique_lock<std::mutex> lck(hashMtx); // we safely collect the new buffers
            buffersVecCpy = buffersVec;
            buffersVec.clear();
        }
        
        for (Buf<uint8_t>* buffers : buffersVecCpy) { // for each array of buffers
            
            for (uint16_t b = 0; b<mapCount; ++b) { // for each buffer file
                
                Buf<uint8_t>* buffer = &buffers[b];
                bufFile[b].write(reinterpret_cast<const char *>(&buffer->pos), sizeof(uint64_t));
                bufFile[b].write(reinterpret_cast<const char *>(buffer->seq), sizeof(uint8_t) * buffer->pos);
                delete[] buffers[b].seq;
                freed += buffers[b].size * sizeof(uint8_t);
                
            }
            
            delete[] buffers;
            
        }
        
    }
    
    for (uint16_t b = 0; b<mapCount; ++b) // we close all files
        bufFile[b].close();
    
    return true;
    
}

bool DBG::buffersToMaps() {
    
    std::vector<std::function<bool()>> jobs;
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of buf files
        fileSizes.push_back(fileSize(userInput.prefix + "/.buf." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        jobs.push_back([this, i] { return processBuffers(i); });
        
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
    consolidateTmpMaps();
    
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

bool DBG::dumpTmpMap(std::string prefix, uint16_t m) {
    
    uint8_t fileNum = 0;
    
    while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin"))
        ++fileNum;
        
    prefix.append("/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    deleteMap(m);
    
    return true;
    
}

void DBG::consolidateTmpMaps(){ // concurrent merging of the maps that store the same hashes
    
    lg.verbose("Consolidating temporary maps");
    
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".0.tmp.bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        mergeTmpMaps(i);
    
}

bool DBG::mergeTmpMaps(uint16_t m) { // a single job merging maps with the same hashes
    
    std::string prefix = userInput.prefix; // loads the first map
    std::string firstFile = prefix + "/.map." + std::to_string(m) + ".0.tmp.bin";
    
    if (!fileExists(prefix + "/.map." + std::to_string(m) + ".1.tmp.bin")) {
        
        std::rename(firstFile.c_str(), (prefix + "/.map." + std::to_string(m) + ".bin").c_str());
        return true;
        
    }
    
    phmap::BinaryInputArchive ar_in(firstFile.c_str());
    maps[m]->phmap_load(ar_in);
    alloc += mapSize(*maps[m]);
    
    remove(firstFile.c_str());
    
    uint8_t fileNum = 0;
    
    while (fileExists(prefix + "/.map." + std::to_string(m) + "." + std::to_string(++fileNum) + ".tmp.bin")) { // for additional map loads the map and merges it
        
        std::string nextFile = prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) + ".tmp.bin"; // loads the next map
        parallelMap* nextMap = new parallelMap;
        phmap::BinaryInputArchive ar_in(nextFile.c_str());
        nextMap->phmap_load(ar_in);
        uint64_t map_size1 = mapSize(*nextMap);
        alloc += map_size1;
        
        uint64_t map_size2 = mapSize(*maps[m]);
        unionSum(nextMap, maps[m], m); // unionSum operation between the existing map and the next map
        
        alloc += mapSize(*maps[m]) - map_size2;
        remove(nextFile.c_str());
        delete nextMap;
        freed += map_size1;
        
    }
    
    dumpMap(userInput.prefix, m);
    
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

bool DBG::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    deleteMap(m);
    
    return true;
    
}

void DBG::finalize() {
    
    if (userInput.inDBG.size() == 0) {
        
        readingDone = true;
        
        joinThreads();
                
        lg.verbose("Converting buffers to maps");
        
        buffersToMaps();
        
    }

}

void DBG::stats() {
    
    lg.verbose("Computing summary statistics");
    
    std::vector<std::function<bool()>> jobs;
    std::array<uint16_t, 2> mapRange = {0,0};
    
    while (mapRange[1] < mapCount) {
        
        mapRange = computeMapRange(mapRange);
        loadMapRange(mapRange);
        
        for (uint32_t i = mapRange[0]; i < mapRange[1]; ++i)
            jobs.push_back([this, i] { return summary(i); });
        
        threadPool.queueJobs(jobs);
        jobWait(threadPool);
        jobs.clear();
        deleteMapRange(mapRange);
    }
    DBGstats();
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

void DBG::DBGstats() {
    
    uint64_t missing = pow(4,k)-totKmersDistinct;
    
    std::cout<<"DBG Summary statistics:\n"
             <<"Total kmers: "<<totKmers<<"\n"
             <<"Unique kmers: "<<totKmersUnique<<"\n"
             <<"Distinct kmers: "<<totKmersDistinct<<"\n"
             <<"Missing kmers: "<<missing<<"\n"
             <<"Total edges: "<<totEdgeCount<<"\n";
    
}

std::array<uint16_t, 2> DBG::computeMapRange(std::array<uint16_t, 2> mapRange) {
    
    uint64_t max = 0;
    mapRange[0] = mapRange[1];
    
    for (uint16_t m = mapRange[0]; m<mapCount; ++m) {
        
        max += fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".bin");
        if(!memoryOk(max))
            break;
        mapRange[1] = m + 1;
        
    }
    
    return mapRange;
    
}

void DBG::loadMapRange(std::array<uint16_t, 2> mapRange) {
    
    std::vector<std::function<bool()>> jobs;
    
    for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
        jobs.push_back([this, m] { return loadMap(userInput.prefix, m); });
    
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    
}

void DBG::deleteMapRange(std::array<uint16_t, 2> mapRange) {
    
    for(uint16_t m = mapRange[0]; m<mapRange[1]; ++m)
        deleteMap(m);
    
}

void DBG::cleanup() {
    
    if(!(userInput.inDBG.size() == 1) && userInput.outFile.find(".kreeq") == std::string::npos) {
        
        lg.verbose("Deleting tmp files");
        
        std::vector<std::function<bool()>> jobs;
        
        for(uint16_t m = 0; m<mapCount; ++m) // remove tmp files
            jobs.push_back([this, m] { return remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str()); });
        
        threadPool.queueJobs(jobs);
        
        if (userInput.prefix != ".")
            rm_dir(userInput.prefix.c_str());
        
    }
    
    jobWait(threadPool);
    
}

void DBG::load() {
    
    if (userInput.inDBG.size() == 1){
        userInput.prefix = userInput.inDBG[0];
    }else if (userInput.inDBG.size() > 1) {
        fprintf(stderr, "More than one DBG database provided. Merge them first. Exiting.\n");
        exit(EXIT_FAILURE);
    }else{
        fprintf(stderr, "Cannot load DBG input. Exiting.\n");
        exit(EXIT_FAILURE);
    }
    
}

bool DBG::loadMap(std::string prefix, uint16_t m) { // loads a specific map
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    phmap::BinaryInputArchive ar_in(prefix.c_str());
//    allocMemory(fileSize(prefix));
    maps[m]->phmap_load(ar_in);
    alloc += mapSize(*maps[m]);
    
    return true;

}

bool DBG::deleteMap(uint16_t m) {
    
    uint64_t map_size = mapSize(*maps[m]);
    delete maps[m];
    freed += map_size;
    
    maps[m] = new parallelMap;
    alloc += mapSize(*maps[m]);
    
    return true;
    
}

bool DBG::updateMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.map." + std::to_string(m) + ".bin");

    if (fileExists(prefix)) {
        
        parallelMap* nextMap = new parallelMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap->phmap_load(ar_in);
        uint64_t map_size1 = mapSize(*nextMap);
        alloc += map_size1;
        
        uint64_t map_size2 = mapSize(*maps[m]);
        unionSum(maps[m], nextMap, m); // merges the current map and the existing map
        alloc += mapSize(*nextMap) - map_size1;
        delete maps[m];
        freed += map_size2;
        maps[m] = nextMap;
    
    }

    dumpMap(userInput.prefix, m);
    return true;
    
}

void DBG::kunion(){ // concurrent merging of the maps that store the same hashes
    
    std::vector<std::function<bool()>> jobs;
    
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.inDBG[0] + "/.map." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        mergeMaps(i);
    
}

bool DBG::mergeMaps(uint16_t m) { // a single job merging maps with the same hashes
    
    std::string prefix = userInput.inDBG[0]; // loads the first map
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
    phmap::BinaryInputArchive ar_in(prefix.c_str());
    maps[m]->phmap_load(ar_in);
    
    unsigned int numFiles = userInput.inDBG.size();

    for (unsigned int i = 1; i < numFiles; ++i) { // for each kmerdb loads the map and merges it
        
        std::string prefix = userInput.inDBG[i]; // loads the next map
        prefix.append("/.map." + std::to_string(m) + ".bin");
        
        parallelMap* nextMap = new parallelMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap->phmap_load(ar_in);
        
        unionSum(nextMap, maps[m], m); // unionSum operation between the existing map and the next map
        delete nextMap;
        
    }
    
    dumpMap(userInput.prefix, m);
    deleteMap(m);
    
    summary(m);
    
    return true;

}

bool DBG::mergeSubMaps(parallelMap* map1, parallelMap* map2, uint8_t subMapIndex, uint16_t m) {
    
    auto& inner = map1->get_inner(subMapIndex);   // to retrieve the submap at given index
    auto& submap1 = inner.set_;        // can be a set or a map, depending on the type of map1
    auto& inner2 = map2->get_inner(subMapIndex);
    auto& submap2 = inner2.set_;
    parallelMap32& map32 = *maps32[m];
    
    for (auto pair : submap1) { // for each element in map1, find it in map2 and increase its value
        
        auto got = submap2.find(pair.first); // insert or find this kmer in the hash table
        if (got == submap2.end()){
            submap2.insert(pair);
        }else{

            DBGkmer& dbgkmerMap = got->second;
        
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
            
        };
        
    }
    return true;
}


bool DBG::unionSum(parallelMap* map1, parallelMap* map2, uint16_t m) {
    
    std::vector<std::function<bool()>> jobs;
    
    if (map1->subcnt() != map2->subcnt()) {
        fprintf(stderr, "Maps don't have the same numbers of submaps (%zu != %zu). Terminating.\n", map1->subcnt(), map2->subcnt());
        exit(EXIT_FAILURE);
    }
    
    for(std::size_t subMapIndex = 0; subMapIndex < map1->subcnt(); ++subMapIndex)
        jobs.push_back([this, map1, map2, subMapIndex, m] { return mergeSubMaps(map1, map2, subMapIndex, m); });
    
    threadPool.queueJobs(jobs);
    jobWait(threadPool);
    
    return true;
    
}
