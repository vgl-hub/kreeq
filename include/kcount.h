#ifndef KCOUNT_H
#define KCOUNT_H

class Kcount {

    unsigned short int k;
    
    unsigned long long int totKmers = 0;
    
    unsigned long long int totKmersUnique = 0;
    
    const unsigned int mapCount = pow(4,k/4);

    std::vector<unsigned long long int>* buff = new std::vector<unsigned long long int>[mapCount];
    
    phmap::flat_hash_map<unsigned long long int, unsigned long long int>* map = new phmap::flat_hash_map<unsigned long long int, unsigned long long int>[mapCount];
    
    const unsigned char ctoi[256] = {
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
    
    Kcount(std::vector<InSegment*>* segments, unsigned short int k) : k(k) {
        
        count(segments);
        
    };
    
    void count(std::vector<InSegment*>* segments);
    
    inline size_t hash(unsigned short int * string);
    
    bool countBuff(std::vector<unsigned long long int>& buff, phmap::flat_hash_map<unsigned long long int, unsigned long long int>& map);
    
    bool countUnique(phmap::flat_hash_map<unsigned long long int, unsigned long long int>& map);
    
    ~Kcount(){
        
        delete[] map;
        delete[] buff;
        
    }

};

#endif /* KCOUNT_H */
