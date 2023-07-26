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
    
    return get_mem_inuse(3) < (userInput.maxMem == 0 ? get_mem_total(3) * 0.4 : userInput.maxMem);
    
}

bool DBG::memoryOk(int64_t delta) {
    
    return get_mem_inuse(3) + convert_memory(delta, 3) < (userInput.maxMem == 0 ? get_mem_total(3) * 0.4 : userInput.maxMem);
    
}

bool DBG::traverseInReads(std::string* readBatch) { // specialized for string objects
    
    {
        std::lock_guard<std::mutex> lck(readMtx);
        readBatches.push(readBatch);
        alloc += readBatch->size() * sizeof(char);
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
    
    int16_t threadN = std::thread::hardware_concurrency() - 2; // the master thread and the writing thread will continue to run so -2
    
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
    
    while (true) {
            
        {
            
            std::unique_lock<std::mutex> lck(readMtx);
            
            if (readingDone && readBatches.size() == 0)
                return true;

            if (readBatches.size() == 0)
                continue;
            
            readBatch = readBatches.front();
            readBatches.pop();
            
        }

        uint64_t len = readBatch->size();
        
        if (len<k) {
            delete readBatch;
            return true;
        }
        
        Buf<kmer> *buffers = new Buf<kmer>[mapCount];
        unsigned char *first = (unsigned char*) readBatch->c_str();
        uint8_t *str = new uint8_t[len];
        uint8_t e = 0;
        uint64_t key, pos = 0, kcount = len-k+1;
        bool isFw = false;
        Buf<kmer>* buffer;
        
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
            pos = buffer->newPos();
            kmer &khmer = buffer->seq[pos];
            khmer.hash = key;
            
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
        
        std::lock_guard<std::mutex> lck(hashMtx);
        for (uint16_t b = 0; b<mapCount; ++b)
            alloc += buffers[b].size * sizeof(kmer);

        freed += len * sizeof(char);
        buffersVec.push_back(buffers);
        
    }
    
    return true;
    
}

bool DBG::dumpBuffers() {
    
    bool hashing = true;
    std::vector<Buf<kmer>*> buffersVecCpy;
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
        
        for (Buf<kmer>* buffers : buffersVecCpy) { // for each array of buffers
            
            for (uint16_t b = 0; b<mapCount; ++b) { // for each buffer file
                
                Buf<kmer>* buffer = &buffers[b];
                bufFile[b].write(reinterpret_cast<const char *>(&buffer->pos), sizeof(uint64_t));
                bufFile[b].write(reinterpret_cast<const char *>(buffer->seq), sizeof(kmer) * buffer->pos);
                delete[] buffers[b].seq;
                freed += buffers[b].size * sizeof(kmer);
                
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
    
    for (uint16_t m = 0; m<mapCount; ++m) // compute size of buffer files
        fileSizes.push_back(fileSize(userInput.prefix + "/.buf." + std::to_string(m) + ".bin"));
    
    std::vector<uint32_t> idx = sortedIndex(fileSizes, true); // sort by largest
    
    for (uint32_t i : idx)
        jobs.push_back([this, i] { return processBuffers(i); });
        
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
    return true;

}

bool DBG::processBuffers(uint16_t m) {
    
    uint64_t pos = 0;
//    int64_t initial_size = 0, final_size = 0;
    Buf<kmer> *buf;
        
    phmap::flat_hash_map<uint64_t, DBGkmer>& map = *maps[m]; // the map associated to this buffer
    
    std::ifstream bufFile(userInput.prefix + "/.buf." + std::to_string(m) + ".bin", std::ios::in | std::ios::binary);
    
    while(bufFile && !(bufFile.peek() == EOF)) {
        
        bufFile.read(reinterpret_cast<char *>(&pos), sizeof(uint64_t));
        
        buf = new Buf<kmer>(pos);
        buf->pos = pos;
        buf->size = pos;
        alloc += buf->size * sizeof(kmer);
        
        bufFile.read(reinterpret_cast<char *>(buf->seq), sizeof(kmer) * buf->pos);
        
        for (uint64_t c = 0; c<pos; ++c) {
            
            kmer &khmer = buf->seq[c];
            
            DBGkmer &dbgkmer = map[khmer.hash];
            
            for (uint64_t w = 0; w<4; ++w) { // update weights
                
                if (255 - dbgkmer.fw[w] >= khmer.fw[w])
                    dbgkmer.fw[w] += khmer.fw[w];
                if (255 - dbgkmer.bw[w] >= khmer.bw[w])
                    dbgkmer.bw[w] += khmer.bw[w];
            }
            if (dbgkmer.cov < 255)
                ++dbgkmer.cov; // increase kmer coverage
            
        }
        
        delete[] buf->seq;
        freed += buf->size * sizeof(kmer);
        delete buf;
        
    }
    
    bufFile.close();
    
    alloc += mapSize(*maps[m]);
    remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str());
    
    dumpMap(userInput.prefix, m);
    
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

void DBG::summary() {
    
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
        jobs.push_back([this, i] { return histogram(i); });
        
    threadPool.queueJobs(jobs);
    
    jobWait(threadPool);
    
    DBGstats();

}

bool DBG::histogram(uint16_t m) {
    
    uint64_t kmersUnique = 0, kmersDistinct = 0, map_size = 0;
    phmap::flat_hash_map<uint64_t, uint64_t> hist;
    
    loadMap(userInput.prefix, m);
    
    for (auto pair : *maps[m]) {
        
        if (pair.second.cov == 1)
            ++kmersUnique;
        
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
    
    for (auto pair : hist) {
        
        finalHistogram[pair.first] += pair.second;
        totKmers += pair.first * pair.second;
        
    }
    
    return true;
    
}

void DBG::DBGstats() {
    
    uint64_t missing = pow(4,k)-totKmersDistinct;
    
    std::cout<<"DBG Summary statistics:\n"
             <<"Total: "<<totKmers<<"\n"
             <<"Unique: "<<totKmersUnique<<"\n"
             <<"Distinct: "<<totKmersDistinct<<"\n"
             <<"Missing: "<<missing<<"\n";
    
}

void DBG::validateSequences(InSequences &inSequences) {
    
    lg.verbose("Validating sequence");
    
    std::vector<InSegment*> *segments = inSequences.getInSegments();
    
    std::array<uint16_t, 2> mapRange = {0,0};
    
    while (mapRange[1] < mapCount-1) {
        
        uint64_t max = 0;
        
        for (uint16_t m = mapRange[0]; m<mapCount; ++m) {
            
            max += fileSize(userInput.prefix + "/.map." + std::to_string(m) + ".bin");

            if(!memoryOk(max))
                break;
            
            mapRange[1] = m;
            
        }
        
        std::vector<std::function<bool()>> jobs;
        
        for(uint16_t m = mapRange[0]; m<=mapRange[1]; ++m)
            jobs.push_back([this, m] { return loadMap(userInput.prefix, m); });
            
        threadPool.queueJobs(jobs);
        
        jobWait(threadPool);
        
        jobs.clear();
        
        for (InSegment* segment : *segments)
            jobs.push_back([this, segment, mapRange] { return validateSegment(segment, mapRange); });
        
        threadPool.queueJobs(jobs);
        
        jobWait(threadPool);
            
        for(uint16_t m = mapRange[0]; m<=mapRange[1]; ++m) {
            
            uint64_t map_size = mapSize(*maps[m]);
            delete maps[m];
            maps[m] = new phmap::flat_hash_map<uint64_t, DBGkmer>;
            freed += map_size;
            
        }
        
        mapRange[0] = mapRange[1] + 1;
        
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
        
        i = key % mapCount;
        
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

void DBG::cleanup() {
    
    if(!(userInput.inDBG.size() == 1) && userInput.outFile == "") {
        
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
    
    alloc += mapSize(*maps[m]);
    
    return true;

}

bool DBG::updateMap(std::string prefix, uint16_t m) {
    
    uint64_t map_size = mapSize(*maps[m]);
    
    prefix.append("/.map." + std::to_string(m) + ".bin");
    
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
    
    histogram(m);
    
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
            
            break;
            
        }
            
    }
    
}
