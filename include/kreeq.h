#ifndef KREEQ_H
#define KREEQ_H

struct kmer {
    
    uint64_t hash = 0;
    uint8_t fw[4] = {0}, bw[4] = {0};
    
};

struct DBGkmer {
    
    uint8_t fw[4] = {0}, bw[4] = {0}, cov = 0;
    
};

class DBG : public Kmap<UserInputKreeq, DBGkmer, kmer> {
    
    std::atomic<uint64_t> totMissingKmers{0}, totKcount{0}, totEdgeMissingKmers{0};
    std::vector<uint32_t> dependencies;
    bool readingDone = false, tmp = false;
    
    UserInputKreeq& userInput;
    
    std::vector<Buf<kmer>*> buffers;

public:
    
    DBG(UserInputKreeq& userInput) : Kmap{userInput.kmerLen} , userInput(userInput) {
        
        lg.verbose("Deleting any tmp file");
        for(uint16_t m = 0; m<mapCount; ++m) // remove tmp files
            threadPool.queueJob([=]{ return remove(("./.kmap." + std::to_string(m) + ".bin").c_str()); });
        
        jobWait(threadPool);
        
        initHashing();
        
    };
    
    std::vector<Log> logs;
    
    bool memoryOk();
    
    bool memoryOk(int64_t delta);
    
    void initHashing();
    
    bool traverseInReads(std::string *readBatch);
    
    bool hashSequences(std::string* readBatch);
    
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
