#ifndef KREEQ_H
#define KREEQ_H

#include <future>

struct kmer {
    
    uint64_t hash = 0;
    bool fw[4] = {0}, bw[4] = {0};
    
};

struct DBGkmer {
    
    uint8_t fw[4] = {0}, bw[4] = {0}, cov = 0;
    
};

class DBG : public Kmap<UserInputKreeq, DBGkmer, kmer> {
    
    std::atomic<uint64_t> totMissingKmers{0}, totKcount{0}, totEdgeMissingKmers{0}, buffers{0};
    std::vector<uint32_t> dependencies;
    bool tmp = false;
    std::atomic<bool> readingDone{false}, dumpMaps{false};
    std::vector<std::thread> threads;
    std::vector<std::future<bool>> futures;
    uint8_t hashThreads = 3;
    std::mutex hashMtx;
    
    UserInputKreeq& userInput;
    
    std::queue<std::string*> readBatches;
    std::vector<uint32_t> buffersDone;
    std::vector<bool> buffingDone = std::vector<bool>(hashThreads, false);

public:
    
    DBG(UserInputKreeq& userInput) : Kmap{userInput.kmerLen} , userInput(userInput) {
        
        lg.verbose("Deleting any tmp file");
        remove((userInput.prefix + "/.buffer.bin").c_str());
        for(uint16_t m = 0; m<mapCount; ++m) // remove tmp maps
            threadPool.queueJob([=]{ return remove((userInput.prefix + "/.kmap." + std::to_string(m) + ".bin").c_str()); });
        
        jobWait(threadPool);
        
        initHashing();
        
    };
    
    std::vector<Log> logs;
    
    bool memoryOk();
    
    bool memoryOk(int64_t delta);
    
    void initHashing();
    
    bool traverseInReads(std::string *readBatch);
    
    bool hashSequences(uint8_t t);
    
    bool processBuffers(std::array<uint16_t, 2> mapRange);
    
    void cleanup();
    
    bool joinBuff(uint16_t m);
    
    void summary();
    
    bool histogram(uint16_t m);
    
    void validateSequences(InSequences &inSequences);
    
    bool validateSegments(std::array<uint16_t, 2> mapRange, std::vector<InSegment*> *segments);
    
    bool validateSegment(InSegment *segment, std::array<uint16_t, 2> mapRange);
    
    void consolidate();
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);
    
    void load();
    
    void report();
    
    bool updateMap(std::string prefix, uint16_t m);
    
    bool unionSum(phmap::flat_hash_map<uint64_t, DBGkmer>& map1, phmap::flat_hash_map<uint64_t, DBGkmer>& map2);
    
};

#endif /* KREEQ_H */
