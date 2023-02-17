#ifndef KREEQ_H
#define KREEQ_H

struct Gkmer {
    
    uint64_t hash = 0;
    unsigned char* base = NULL;
    unsigned int sUId = 0, cov = 0;
    
};

class Kpos : public Kmap<UserInputKreeq, std::vector<Gkmer*>, Gkmer> {

    using Kmap<UserInputKreeq, std::vector<Gkmer*>, Gkmer>::Kmap;
    
public:
   
    bool traverseInReads(Sequences* readBatch);
    
    void hashSegments();
    
    void hashSequences(Sequences* readBatch);
    
    void index();
    
    bool joinBuff(uint16_t m);
    
    bool histogram(phmap::flat_hash_map<uint64_t, std::vector<Gkmer*>>& map);
    
    void report(UserInputKreeq userInput);
    
};

struct DBGkmer {
    
    uint64_t hash = 0;
    bool fw[4] = {false}, bw[4] = {false};
    unsigned int cov = 0;
    
};

class DBG : public Kmap<UserInputKreeq, DBGkmer, DBGkmer> {

    using Kmap<UserInputKreeq, DBGkmer, DBGkmer>::Kmap;
    
    uint64_t totErrorKmers = 0, totKcount = 0;

public:
    
    std::vector<Log> logs;
   
    bool traverseInReads(Sequences* readBatch);
    
    void hashSequences(Sequences* readBatch);
    
    void build();
    
    bool joinBuff(uint16_t m);
    
    bool histogram(uint16_t m);
    
    void validateSequences(InSequences &inSequences);
    
    bool validateSegment(InSegment* segment);
    
};

#endif /* KREEQ_H */
