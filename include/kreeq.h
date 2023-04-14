#ifndef KREEQ_H
#define KREEQ_H

struct DBGkmer {
    
    uint64_t hash = 0;
    bool fw[4] = {false}, bw[4] = {false};
    unsigned int cov = 0;
    
};

class DBG : public Kmap<UserInputKreeq, DBGkmer, DBGkmer> {
    
    uint64_t totErrorKmers = 0, totKcount = 0;

public:
    
    DBG(uint8_t k) : Kmap{k} {};
    
    std::vector<Log> logs;
    
    void hashSequences(std::string* readBatch);
    
    void finalize();
    
    bool joinBuff(uint16_t m);
    
    bool histogram(uint16_t m);
    
    void validateSequences(InSequences &inSequences);
    
    bool validateSegment(InSegment* segment);
    
    bool countBuffs(uint16_t m);
    
    bool countBuff(Buf<DBGkmer>* thisBuf, uint16_t m);
    
    void consolidate();
    
    bool traverseInReads(std::string* readBatch);
    
};

#endif /* KREEQ_H */
