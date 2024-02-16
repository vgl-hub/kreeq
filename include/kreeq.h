#ifndef KREEQ_H
#define KREEQ_H

#include <future>

struct edgeBit {
    
    uint8_t edges = 0;
    
    void assign(const uint8_t edge) {
        edges |= 1 << (7 - edge);
    }
    
    bool read(const uint8_t edge) {
        return edges & 1 << (7 - edge);
    }
    
};

struct DBGkmer {
    
    uint8_t fw[4] = {0}, bw[4] = {0}, cov = 0;
    
};

using parallelMap = phmap::parallel_flat_hash_map<uint64_t, DBGkmer,
                                          std::hash<uint64_t>,
                                          std::equal_to<uint64_t>,
                                          std::allocator<std::pair<const uint64_t, DBGkmer>>,
                                          8,
                                          phmap::NullMutex>;


class DBG : public Kmap<UserInputKreeq, DBGkmer, uint8_t> {
    
    std::atomic<uint64_t> totMissingKmers{0}, totKcount{0}, totEdgeMissingKmers{0}, buffers{0};
    std::atomic<bool> readingDone{false};
    std::vector<std::thread> threads;
    std::vector<std::future<bool>> futures;
    std::mutex readMtx, hashMtx;
    std::chrono::high_resolution_clock::time_point past;
    
    UserInputKreeq &userInput;
    InSequencesDBG *genome;
    
    std::queue<std::string*> readBatches;
    
    uint64_t totEdgeCount = 0;

public:
    
    DBG(UserInputKreeq& userInput) : Kmap{userInput.kmerLen} , userInput(userInput) {
        
        lg.verbose("Deleting any tmp file");
        for(uint16_t m = 0; m<mapCount; ++m) {// remove tmp buffers and maps if any
            threadPool.queueJob([=]{ return remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str()); });
            threadPool.queueJob([=]{ return remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str()); });
            uint8_t fileNum = 0;
            while (fileExists(userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum++) +  ".tmp.bin"))
                threadPool.queueJob([=]{ return remove((userInput.prefix + "/.map." + std::to_string(m) + "." + std::to_string(fileNum) +  ".tmp.bin").c_str()); });
            remove((userInput.prefix + "/.index").c_str());
        }
            
        jobWait(threadPool);
        
        if (userInput.inDBG.size() == 0) // start parallel hashing
            initHashing();
        
    };
    
    std::vector<Log> logs;
    
    void status();
    
    void joinThreads();
    
    bool memoryOk();
    
    bool memoryOk(int64_t delta);
    
    bool traverseInReads(std::string *readBatch);
    
    void consolidate();
    
    void initHashing();
    
    bool hashSequences();
    
    bool dumpBuffers();
    
    bool buffersToMaps();
    
    bool processBuffers(uint16_t m);
    
    void cleanup();
    
    bool joinBuff(uint16_t m);
    
    void finalize();
    
    bool summary(uint16_t m);
    
    void DBGstats();
    
    void loadGenome(InSequencesDBG *genome);
    
    void validateSequences();
    
    bool evaluateSegment(uint32_t i, std::array<uint16_t, 2> mapRange);
    
    bool validateSegment(uint32_t i);
    
    bool dumpTmpMap(std::string prefix, uint16_t m);
    
    void consolidateTmpMaps();
    
    bool mergeTmpMaps(uint16_t m);
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);
    
    bool deleteMap(uint16_t m);
    
    void load();
    
    bool updateMap(std::string prefix, uint16_t m);
    
    void kunion();
    
    bool mergeMaps(uint16_t m);
    
    bool mergeSubMaps(parallelMap* map1, parallelMap* map2, uint8_t subMapIndex);
    
    bool unionSum(parallelMap* map1, parallelMap* map2);
    
    void report();
    
    void printTable(std::string ext);
    
    void printTableCompressed();
    
    void printTableCompressedBinary();
    
    std::array<uint16_t, 2> computeMapRange(std::array<uint16_t, 2> mapRange);
    
    void loadMapRange(std::array<uint16_t, 2> mapRange);
    
    void deleteMapRange(std::array<uint16_t, 2> mapRange);
    
};

#endif /* KREEQ_H */
