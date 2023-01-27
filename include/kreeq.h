#ifndef KREEQ_H
#define KREEQ_H

struct Gkmer {
    
    uint64_t hash = 0;
    unsigned char* base = NULL;
    unsigned int sUId = 0, cov = 0;
    
};

class Kpos : public Kmap<UserInputKreeq, std::vector<Gkmer*>, Gkmer> {

    using Kmap<UserInputKreeq, std::vector<Gkmer*>, Gkmer>::Kmap;
    
    friend class InSequences;

public:
   
    bool traverseInReads(Sequences* readBatch);
    
    void hashSegments();
    
    void index();
    
    bool joinBuff(uint16_t m);
    
    bool histogram(phmap::flat_hash_map<uint64_t, std::vector<Gkmer*>>& map);
    
    void report(UserInputKreeq userInput);
    
};

#endif /* KREEQ_H */
