#ifndef KREEQ_H
#define KREEQ_H

class pos {
    
    char* pos = NULL;
    unsigned int cov = 0;
    
};

class Kpos : public Kmap<UserInputKreeq, std::vector<pos>> {
    
    using Kmap<UserInputKreeq, std::vector<pos>>::Kmap;

public:
    
    bool validate(UserInputKreeq userInput);
    
    bool joinBuff(uint16_t m);
    
    bool stats(phmap::flat_hash_map<uint64_t, std::vector<pos>>& map);
    
};

#endif /* KREEQ_H */
