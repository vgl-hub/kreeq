#ifndef KREEQ_H
#define KREEQ_H

class pos {
    
    char* pos = NULL;
    unsigned int cov = 0;
    
};

class Kpos : public Kmap<UserInputKreeq, std::vector<pos>> {
    
    InSequences inSequences;
    
    using Kmap<UserInputKreeq, std::vector<pos>>::Kmap;
    
    friend class InSequences;

public:
    
    bool convert(UserInputKreeq userInput);
   
    bool traverseInReads(Sequences* readBatch);
    
    void hashSequences(Sequences* readBatch);
    
    void validate(UserInputKreeq userInput);
    
    bool joinBuff(uint16_t m);
    
    bool stats(phmap::flat_hash_map<uint64_t, std::vector<pos>>& map);
    
    void report(UserInputKreeq userInput);
    
};

#endif /* KREEQ_H */
