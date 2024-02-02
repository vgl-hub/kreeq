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

#include "input.h"
#include "kmer.h"
#include "kreeq.h"

uint64_t mapSize(phmap::flat_hash_map<uint64_t, DBGkmer>& m) {
    
    return (m.size() * (sizeof(DBGkmer) + sizeof(void*)) + // data list
     m.bucket_count() * (sizeof(void*) + sizeof(uint64_t))) // bucket index
    * 1.3; // estimated allocation overheads
    
}

double errorRate(uint64_t missingKmers, uint64_t totalKmers, uint8_t k){ // estimate QV from missing kmers
    
    return 1 - pow(1 - (double) missingKmers/totalKmers, (double) 1/k);
    
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
    
    allocMemory(readBatch->size() * sizeof(char));
    
    {
        std::lock_guard<std::mutex> lck(readMtx);
        readBatches.push(readBatch);
    }
    
    return true;
    
}

void DBG::consolidate() {
    
    status();
    
    while (!memoryOk()){status();}

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
            
            readBatch = readBatches.front();
            readBatches.pop();
            len = readBatch->size();
            
        }
        
        allocMemory(len * sizeof(uint8_t)); // this is for the string we allocate next
        
        if (len<k) {
            delete readBatch;
            continue;
        }
        
        for (uint16_t b = 0; b<mapCount; ++b)
            allocMemory(pow(2,8) * sizeof(uint8_t));
        Buf<uint8_t> *buffers = new Buf<uint8_t>[mapCount];
        unsigned char *first = (unsigned char*) readBatch->c_str();
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
    
    for(uint16_t b = 0; b<mapCount; ++b) {
        std::string fl = userInput.prefix + "/.buf." + std::to_string(b) + ".bin";
        allocMemory(fileSize(fl));
        threadPool.queueJob([this, b] { return processBuffers(b); });
    }
    
    jobWait(threadPool);
    
    return true;

}

bool DBG::processBuffers(uint16_t m) {
    
    uint64_t pos = 0, hash;
    Buf<uint8_t> *buf;
    edgeBit edges;
    
    std::string fl = userInput.prefix + "/.buf." + std::to_string(m) + ".bin";
    uint64_t flSize = fileSize(fl);
    std::ifstream bufFile(fl, std::ios::in | std::ios::binary);
    phmap::flat_hash_map<uint64_t, DBGkmer>& map = *maps[m]; // the map associated to this buffer
    map.reserve(flSize / 17); // 8 + 8 + 1
    
    while(bufFile && !(bufFile.peek() == EOF)) {
        
        bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
        
        allocMemory(pos * sizeof(uint8_t));
        buf = new Buf<uint8_t>(pos);
        
        buf->pos = pos;
        buf->size = pos;
        
        bufFile.read(reinterpret_cast<char *>(buf->seq), sizeof(uint8_t) * buf->pos);
        
        for (uint64_t c = 0; c<pos; c+=9) {
            
            memcpy(&hash, &buf->seq[c], 8);
            memcpy(&edges, &buf->seq[c+8], 1);
            
            DBGkmer &dbgkmer = map[hash];
            
            for (uint64_t w = 0; w<4; ++w) { // update weights
                
                if (255 - dbgkmer.fw[w] >= edges.read(w))
                    dbgkmer.fw[w] += edges.read(w);
                if (255 - dbgkmer.bw[w] >= edges.read(4+w))
                    dbgkmer.bw[w] += edges.read(4+w);
            }
            if (dbgkmer.cov < 255)
                ++dbgkmer.cov; // increase kmer coverage
            
        }
        
        delete[] buf->seq;
        freed += buf->size * sizeof(uint8_t);
        delete buf;
        
    }
    
    alloc += mapSize(*maps[m]) - flSize;
    
    dumpMap(userInput.prefix, m); // if it does, dump map
    
    bufFile.close();
    
    remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str());
    
    return true;
    
}

bool DBG::dumpMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
    phmap::BinaryOutputArchive ar_out(prefix.c_str());
    maps[m]->phmap_dump(ar_out);
    
    uint64_t map_size = mapSize(*maps[m]);
    delete maps[m];
    freed += map_size;
    
    maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
    
    return true;
    
}

void DBG::finalize() {
    
    if (userInput.inDBG.size() == 0) {
        
        readingDone = true;
        
        joinThreads();
                
        lg.verbose("Loading buffers in maps");
        
        buffersToMaps();
        
    }
    
    lg.verbose("Computing summary statistics");
    
    std::vector<std::function<bool()>> jobs;
    
    std::vector<uint64_t> fileSizes;
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of map files
        fileSizes.push_back(fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for (uint32_t i : idx)
        jobs.push_back([this, i] { return summary(i); });
        
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
    DBGstats();

}

bool DBG::summary(uint16_t m) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0, map_size = 0, edgeCount = 0;
    phmap::flat_hash_map<uint64_t, uint64_t> hist;
    
    loadMap(userInput.prefix, m);
    
    for (auto pair : *maps[m]) {
        
        if (pair.second.cov == 1)
            ++kmersUnique;
        
        for (uint64_t w = 0; w<4; ++w) // update weights
            edgeCount += pair.second.fw[w] > 0 ? 1 : 0 + pair.second.bw[w] > 0 ? 1 : 0;
        
        ++kmersDistinct;
        ++hist[pair.second.cov];
        
    }
    
    map_size = mapSize(*maps[m]);
    delete maps[m];
    maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
    freed += map_size;
 
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

void DBG::loadGenome(InSequencesDBG *genome) {
    
    this->genome = genome;
}

void DBG::validateSequences() {
    
    std::vector<std::function<bool()>> jobs;
    
    lg.verbose("Validating sequence");
    
    genome->generateValidationVector();
    
    std::vector<InSegment*> *segments = genome->getInSegments();
    
    std::array<uint16_t, 2> mapRange = {0,0};
    
    while (mapRange[1] < mapCount-1) {
        
        uint64_t max = 0;
        
        for (uint16_t m = mapRange[0]; m<mapCount; ++m) {
            
            max += fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".bin");

            if(!memoryOk(max))
                break;
            
            mapRange[1] = m;
            
        }
        
        for(uint16_t m = mapRange[0]; m<=mapRange[1]; ++m)
            jobs.push_back([this, m] { return loadMap(userInput.prefix, m); });
            
        threadPool.queueJobs(jobs);
        
        jobWait(threadPool);
        
        jobs.clear();
        
        for (uint32_t s = 0; s < segments->size(); ++s)
            jobs.push_back([this, s, mapRange] { return evaluateSegment(s, mapRange); });
        
        threadPool.queueJobs(jobs);
        
        jobWait(threadPool);
        
        jobs.clear();
            
        for(uint16_t m = mapRange[0]; m<=mapRange[1]; ++m) {
            
            uint64_t map_size = mapSize(*maps[m]);
            delete maps[m];
            maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
            freed += map_size;
            
        }
        
        mapRange[0] = mapRange[1] + 1;
        
    }
    
    jobs.clear();
    
    for (uint32_t s = 0; s < segments->size(); ++s)
        jobs.push_back([this, s] { return validateSegment(s); });
    
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
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

bool DBG::evaluateSegment(uint32_t s, std::array<uint16_t, 2> mapRange) {
    
    std::vector<uint64_t> missingKmers;
    std::vector<uint64_t> edgeMissingKmers;
    
    std::vector<InSegment*> *segments = genome->getInSegments();
    std::vector<DBGbase*> *dbgbases = genome->getInSegmentsDBG();
    
    InSegment *segment = (*segments)[s];
    DBGbase *DBGsequence = (*dbgbases)[s];
    
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
        
        i = key % mapCount;
        
//        std::cout<<"\n"<<itoc[*(str+c)]<<"\t"<<c<<"\t"<<isFw<<std::endl;
        
        if (i >= mapRange[0] && i <= mapRange[1]) {
            
            map = maps[i];
            
            auto it = map->find(key);
            
            DBGkmer khmer;
            
            if (it != map->end()) {
                khmer = it->second;
                DBGsequence[c].cov = khmer.cov;
                DBGsequence[c].isFw = isFw;
            }
            
            if (DBGsequence[c].cov == 0) // merqury QV
                missingKmers.push_back(c);
            else if (DBGsequence[c].cov < userInput.covCutOff) // merqury QV with cutoff
                missingKmers.push_back(c);
            else { // kreeq QV
                bool noEdgeLeft = false, noEdgeRight = false;
                if (DBGsequence[c].isFw){

                    if (c<kcount-1) {
                        if (khmer.fw[*(str+c+k)] != 0)
                            DBGsequence[c].fw = khmer.fw[*(str+c+k)];
                        else
                            noEdgeRight = true;
                    }
                    if (c>0) {
                        if (khmer.bw[*(str+c-1)] != 0)
                            DBGsequence[c].bw = khmer.bw[*(str+c-1)];
                        else
                            noEdgeLeft = true;
                        // std::cout<<"edge error1"<<std::endl;
                    }
                    
                }else{
                    
                    if (c>0) {
                        if (khmer.fw[3-*(str+c-1)] != 0)
                            DBGsequence[c].fw = khmer.fw[3-*(str+c-1)];
                        else
                            noEdgeLeft = true;
                    }
                    
                    if (c<kcount-1) {
                        if (khmer.bw[3-*(str+c+k)] != 0)
                            DBGsequence[c].bw = khmer.bw[3-*(str+c+k)];
                        else
                            noEdgeRight = true;
                        // std::cout<<"edge error2"<<std::endl;
                    }
                }
                if (noEdgeLeft && noEdgeRight) {
                    edgeMissingKmers.push_back(c);
//                    std::cout<<"edge error"<<std::endl;
                }
            }
            ++kmers;

        }
    }
    
    delete[] str;
    
    totMissingKmers += missingKmers.size();
    totKcount += kmers;
    totEdgeMissingKmers += edgeMissingKmers.size();
    
    return true;
    
}

bool DBG::validateSegment(uint32_t s) {
    
    // internal
    uint32_t a = s;
    s = a;
    
    return true;
    
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
    maps[m]->phmap_load(ar_in);
    
    allocMemory(mapSize(*maps[m]));
    
    return true;

}

bool DBG::updateMap(std::string prefix, uint16_t m) {
    
    prefix.append("/.map." + std::to_string(m) + ".bin");

    if (fileExists(prefix)) {
        
        phmap::flat_hash_map<uint64_t, DBGkmer> *nextMap = new phmap::flat_hash_map<uint64_t, DBGkmer>;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        allocMemory(fileSize(prefix));
        nextMap->phmap_load(ar_in);
        freed += fileSize(prefix);
        alloc += mapSize(*nextMap);
    
        uint64_t map_size1 = mapSize(*nextMap);
        uint64_t map_size2 = mapSize(*maps[m]);
        unionSum(*maps[m], *nextMap); // merges the current map and the existing map
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
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of buffer files
        fileSizes.push_back(fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for(uint32_t i : idx)
        jobs.push_back([this, i] { return mergeMaps(i); });
    
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
    DBGstats();
    
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
        
        phmap::flat_hash_map<uint64_t, DBGkmer> nextMap;
        phmap::BinaryInputArchive ar_in(prefix.c_str());
        nextMap.phmap_load(ar_in);
        
        unionSum(nextMap, *maps[m]); // unionSum operation between the existing map and the next map
        
    }
    
    dumpMap(userInput.prefix, m);
    uint64_t map_size = mapSize(*maps[m]);
    delete maps[m];
    maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
    freed += map_size;
    
    summary(m);
    
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

void DBG::report() { // generates the output from the program
    
    const static phmap::flat_hash_map<std::string,int> string_to_case{ // different outputs available
        {"kreeq",1},
        {"bed",2}
    };
    
    std::string ext = "stdout";
    
    if (userInput.outFile != "")
        ext = getFileExt("." + userInput.outFile);
    
    lg.verbose("Writing ouput: " + userInput.outFile);
    
    std::unique_ptr<std::ostream> ostream; // smart pointer to handle any kind of output stream
    
    switch (string_to_case.count(ext) ? string_to_case.at(ext) : 0) {
        
        default:
        case 1: { // .kreeq
            
            std::ofstream ofs(userInput.outFile + "/.index"); // adding index
            
            ostream = std::make_unique<std::ostream>(ofs.rdbuf());
            
            *ostream<<+k<<"\n"<<mapCount<<std::endl;
            
            ofs.close();
            
            break;
            
        }
            
        case 2: { // .bed
            
            printTable();
            
            break;
            
        }
            
    }
    
}

void DBG::printTable() {
    
    std::ofstream ofs(userInput.outFile);
    
    genome->sortPathsByOriginal();
    
    std::vector<InPath> inPaths = genome->getInPaths();
    std::vector<InSegment*> *inSegments = genome->getInSegments();
    std::vector<InGap> *inGaps = genome->getInGaps();
    std::vector<DBGbase*> *dbgbases = genome->getInSegmentsDBG();
    
    for (InPath& path : inPaths) {
        
        unsigned int cUId = 0, gapLen = 0, sIdx = 0;
        
        std::vector<PathComponent> pathComponents = path.getComponents();
        
        uint64_t absPos = 0;
        
        for (std::vector<PathComponent>::iterator component = pathComponents.begin(); component != pathComponents.end(); component++) {
            
            cUId = component->id;
            
            if (component->type == SEGMENT) {
                
                auto inSegment = find_if(inSegments->begin(), inSegments->end(), [cUId](InSegment* obj) {return obj->getuId() == cUId;}); // given a node Uid, find it
                
                if (inSegment != inSegments->end()) {sIdx = std::distance(inSegments->begin(), inSegment);} // gives us the segment index
                
                DBGbase *dbgbase = (*dbgbases)[sIdx];
                
                std::vector<uint8_t> kmerCov(k-1,0);
                std::vector<uint8_t> edgeCovFw(k-1,0);
                std::vector<uint8_t> edgeCovBw(k-1,0);
                
                if (component->orientation == '+') {
                    
                    for (uint64_t i = 0; i < (*inSegment)->getSegmentLen(); ++i) {
                        
                        ofs<<path.getHeader()
                           <<"\t"<<absPos<<"\t";
                        
                        kmerCov.push_back(dbgbase[i].cov);
                        
                        for(uint8_t c = 0; c<k; ++c){
                            ofs<<std::to_string(kmerCov[c]);
                            if (c < k - 1)
                                ofs<<":";
                        }
                        
                        ofs<<",";
                        
                        edgeCovFw.push_back(dbgbase[i].isFw ? dbgbase[i].fw : dbgbase[i].bw);
                        
                        for(uint8_t c = 0; c<k; ++c){
                            ofs<<std::to_string(edgeCovFw[c]);
                            if (c < k - 1)
                                ofs<<":";
                        }
                        
                        ofs<<",";
                        
                        edgeCovBw.push_back(dbgbase[i].isFw ? dbgbase[i].bw : dbgbase[i].fw);
                        
                        for(uint8_t c = 0; c<k; ++c){
                            ofs<<std::to_string(edgeCovBw[c]);
                            if (c < k - 1)
                                ofs<<":";
                        }
                        
                        ofs<<"\n";
                        
                        kmerCov.erase(kmerCov.begin());
                        edgeCovFw.erase(edgeCovFw.begin());
                        edgeCovBw.erase(edgeCovBw.begin());
                        
                        ++absPos;
                        
                    }
                    
                }else{
                    
                    // GFA not handled yet
                    
                }
                
                kmerCov.clear();
                
            }else if (component->type == GAP){
                
                auto inGap = find_if(inGaps->begin(), inGaps->end(), [cUId](InGap& obj) {return obj.getuId() == cUId;}); // given a node Uid, find it
                
                gapLen += inGap->getDist(component->start - component->end);
                
                absPos += gapLen;
                
            }else{} // need to handle edges, cigars etc
            
        }
        
//        for (uint64_t p = 0; p<kcount; ++p) {
//
//            for (uint8_t c = e; c<k; ++c) { // generate k bases if e=0 or the next if e=k-1
//
//                str[p+c] = ctoi[*(first+p+c)]; // convert the next base
//                if (str[p+c] > 3) { // if non-canonical base is found
//                    p = p+c; // move position
//                    e = 0; // reset base counter
//                    break;
//                }
//
//                e = k-1;
//
//            }
//
//        }
        
    }
    
    ofs.close();
    
}
