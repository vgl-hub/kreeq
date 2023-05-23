#ifndef KREEQ_H
#define KREEQ_H

struct DBGkmer {
    
    uint64_t hash = 0;
    uint8_t cov = 0, fw[4] = {0}, bw[4] = {0};
    
};

class DBG : public Kmap<UserInputKreeq, DBGkmer, DBGkmer> {
    
    uint64_t totMissingKmers = 0, totKcount = 0;
    std::vector<uint32_t> dependencies;
    bool tmp = false;
    
    UserInputKreeq& userInput;

public:
    
    DBG(UserInputKreeq& userInput) : Kmap{userInput.kmerLen} , userInput(userInput) {};
    
    std::vector<Log> logs;
    
    bool hashSequences(std::string *readBatch);
    
    void finalize();
    
    bool joinBuff(uint16_t m);
    
    bool histogram(uint16_t m);
    
    void validateSequences(InSequences &inSequences);
    
    bool validateSegment(InSegment *segment);
    
    bool countBuffs(uint16_t m);
    
    bool countBuff(Buf<DBGkmer> *thisBuf, uint16_t m);
    
    void consolidate();
    
    bool traverseInReads(std::string *readBatch);
    
    bool dumpMap(std::string prefix, uint16_t m);
    
    bool loadMap(std::string prefix, uint16_t m);
    
    void load();
    
    void report();
    
    void updateDBG();
    
    bool updateMap(std::string prefix, uint16_t m);
    
    bool unionSum(phmap::flat_hash_map<uint64_t, DBGkmer>& map1, phmap::flat_hash_map<uint64_t, DBGkmer>& map2);
    
};

#endif /* KREEQ_H */
