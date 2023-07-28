#ifndef KREEQ_H
#define KREEQ_H

#include <future>

struct kmer {
    
    uint64_t hash = 0;
    bool fw[4] = {0}, bw[4] = {0};
    
};

struct edgeBit {
    
    uint8_t edges = 0;
    
    void assign(const uint8_t edge) {
        edges |= 1 << (7 - edge);
    }
    
    bool read(const uint8_t edge) {
        return edges & 1 << (7 - edge);
    }
    
};

struct kmer2 {
    
    uint64_t hash = 0;
    edgeBit edges;
    
};

struct DBGkmer {
    
    uint8_t fw[4] = {0}, bw[4] = {0}, cov = 0;
    
};

class DBG : public Kmap<UserInputKreeq, DBGkmer, kmer> {
    
    std::atomic<uint64_t> totMissingKmers{0}, totKcount{0}, totEdgeMissingKmers{0}, buffers{0};
    std::atomic<bool> readingDone{false};
    std::vector<std::thread> threads;
    std::vector<std::future<bool>> futures;
    std::mutex readMtx, hashMtx;
    std::chrono::high_resolution_clock::time_point past;
    
    UserInputKreeq& userInput;
    
    std::queue<std::string*> readBatches;
    std::vector<Buf<uint8_t>*> buffersVec;

public:
    
    DBG(UserInputKreeq& userInput) : Kmap{userInput.kmerLen} , userInput(userInput) {
        
        lg.verbose("Deleting any tmp file");
        for(uint16_t m = 0; m<mapCount; ++m) {// remove tmp buffers and maps if any
            threadPool.queueJob([=]{ return remove((userInput.prefix + "/.map." + std::to_string(m) + ".bin").c_str()); });
            threadPool.queueJob([=]{ return remove((userInput.prefix + "/.buf." + std::to_string(m) + ".bin").c_str()); });
        }
            
        jobWait(threadPool);
        
        if (userInput.inDBG.size() == 0)
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
    
    void summary();
    
    bool histogram(uint16_t m);
    
    void DBGstats();
    
    void validateSequences(InSequences &inSequences);
    
    bool validateSegments(std::array<uint16_t, 2> mapRange, std::vector<InSegment*> *segments);
    
    bool validateSegment(InSegment *segment, std::array<uint16_t, 2> mapRange);
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);
    
    void load();
    
    void report();
    
    bool updateMap(std::string prefix, uint16_t m);
    
    void kunion();
    
    bool mergeMaps(uint16_t m);
    
    bool unionSum(phmap::flat_hash_map<uint64_t, DBGkmer>& map1, phmap::flat_hash_map<uint64_t, DBGkmer>& map2);
    
};

#endif /* KREEQ_H */
