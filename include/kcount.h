#ifndef KCOUNT_H
#define KCOUNT_H

struct buf64 {
    uint64_t pos = 0, size = 100000;
    uint64_t *seq = new uint64_t[size];
};

class Kcount {

    uint8_t k;
    
    uint64_t totKmers = 0;
    
    uint64_t totKmersUnique = 0;
    
    const uint64_t mapCount = k < 28 ? pow(4,k/4) : pow(4,6);

    buf64* buf = new buf64[mapCount];
    
    phmap::flat_hash_map<uint64_t, uint64_t>* map = new phmap::flat_hash_map<uint64_t, uint64_t>[mapCount];
    
    const uint8_t ctoi[256] = {
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 0, 4, 1, 4, 4, 4, 2, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 3, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4,
          4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4, 4
    };
    
public:
    
    Kcount(std::vector<InSegment*>* segments, uint8_t k) : k(k) {
        
        count(segments);
        
    };
    
    void count(std::vector<InSegment*>* segments);
    
    inline uint64_t hash(uint8_t* string);
    
    bool countBuff(buf64* buf, phmap::flat_hash_map<uint64_t, uint64_t>& map);
    
    bool countUnique(phmap::flat_hash_map<uint64_t, uint64_t>& map);
    
    void resizeBuff(buf64* buff);
    
    ~Kcount(){
        
        delete[] map;
        
        for(uint16_t i = 0; i<mapCount; i++)
            delete[] buf[i].seq;
        
        delete[] buf;
        
    }

};

#endif /* KCOUNT_H */
