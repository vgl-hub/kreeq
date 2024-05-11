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

#define LARGEST 4294967295 // 2^32-1
struct DBGkmer32 {
    uint32_t fw[4] = {0}, bw[4] = {0}, cov = 0;
    
    DBGkmer32() {}
    
    DBGkmer32(const DBGkmer& dbgkmer) {
        std::copy(std::begin(dbgkmer.fw), std::end(dbgkmer.fw), std::begin(fw));
        std::copy(std::begin(dbgkmer.bw), std::end(dbgkmer.bw), std::begin(bw));
        cov = dbgkmer.cov;
    }
    
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
    
    using parallelMap32 = phmap::parallel_flat_hash_map<uint64_t, DBGkmer32,
                                              std::hash<uint64_t>,
                                              std::equal_to<uint64_t>,
                                              std::allocator<std::pair<const uint64_t, DBGkmer32>>,
                                              8,
                                              phmap::NullMutex>;
    
    std::vector<parallelMap32*> maps32;

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
        
        for(uint16_t m = 0; m<mapCount; ++m)
            maps32.push_back(new parallelMap32);
        
    };
    
    ~DBG(){ // always need to call the destructor and delete for any object called with new to avoid memory leaks
        for (parallelMap32* map : maps32)
            delete map;
    }
    
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
    
    void stats();
    
    bool summary(uint16_t m);
    
    void DBGstats();
    
    void loadGenome(InSequencesDBG *genome);
    
    void correctSequences();
    
    bool evaluateSegment(uint32_t i, std::array<uint16_t, 2> mapRange);
    
    void validateSequences();
    
    bool dumpTmpMap(std::string prefix, uint16_t m);
    
    void consolidateTmpMaps();
    
    bool mergeTmpMaps(uint16_t m);
    
    bool reloadMap32(uint16_t m);
    
    bool dumpHighCopyKmers();
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);
    
    bool loadHighCopyKmers();
    
    bool deleteMap(uint16_t m);
    
    void load();
    
    bool updateMap(std::string prefix, uint16_t m);
    
    void kunion();
    
    bool mergeMaps(uint16_t m);
    
    bool mergeSubMaps(parallelMap* map1, parallelMap* map2, uint8_t subMapIndex, uint16_t m);
    
    bool unionSum(parallelMap* map1, parallelMap* map2, uint16_t m);
    
    void report();
    
    void printTable(std::string ext);
    
    void printTableCompressed();
    
    void writeIndex(std::ofstream &ofs);
    
    void printTableCompressedBinary();
    
    void printGFA();
    
    void printVCF();
    
    bool searchGraph(std::array<uint16_t, 2> mapRange);
    
    std::pair<DBGkmer*,bool> findDBGkmer(uint8_t *origin);
    
    int8_t scorePath(std::string path);
    
    int8_t checkNext(uint8_t *nextKmer, uint8_t *origin);
    
    std::deque<DBGpath> findPaths(uint8_t *origin, uint8_t *target, uint8_t depth, DBGpath currentPath, Log &threadLog);
    
    void printAltPaths(std::vector<std::vector<uint8_t>> altPaths, Log &threadLog);
    
    bool loadAnomalies(InSegment *inSegment, std::vector<uint64_t> &anomalies);
    
    BedCoordinates BEDPathsToSegments();
    
    bool detectAnomalies(InSegment *inSegment, std::vector<uint64_t> &anomalies);
    
    bool DBGtoVariants(InSegment *inSegment);
    
    bool variantsToGFA(InSegment *inSegment, Log &threadLog);
    
    std::array<uint16_t, 2> computeMapRange(std::array<uint16_t, 2> mapRange);
    
    void loadMapRange(std::array<uint16_t, 2> mapRange);
    
    void deleteMapRange(std::array<uint16_t, 2> mapRange);
    
};

#endif /* KREEQ_H */
