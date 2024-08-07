#ifndef KREEQ_H
#define KREEQ_H

#include <future>

struct edgeBit {
    
    uint8_t edges = 0;
    
    void assign(const uint8_t edge) {
        edges |= 1 << (7 - edge);
    }
    
    bool read(const uint8_t edge) {
        return edges & 1 << (7 - edge);
    }
    
};

struct DBGkmer {
    uint8_t fw[4] = {0}, bw[4] = {0}, cov = 0;
    
    uint8_t fwCount() { // count forward edges
        
        uint8_t fwEdges = 0;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (fw[i] != 0)
                ++fwEdges;
        }
        return fwEdges;
    }
    
    uint8_t bwCount() { // count backward edges
        
        uint8_t bwEdges = 0;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (bw[i] != 0)
                ++bwEdges;
        }
        return bwEdges;
    }
    
    std::vector<uint8_t> fwEdgeIndexes() { // count forward edges
        
        std::vector<uint8_t> indexes;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (fw[i] != 0)
                indexes.push_back(i);
        }
        return indexes;
    }
    
    std::vector<uint8_t> bwEdgeIndexes() { // count forward edges
        
        std::vector<uint8_t> indexes;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (bw[i] != 0)
                indexes.push_back(i);
        }
        return indexes;
    }
};

#define LARGEST 4294967295 // 2^32-1
struct DBGkmer32 {
    uint32_t fw[4] = {0}, bw[4] = {0}, cov = 0;
    
    DBGkmer32() {}
    
    DBGkmer32(const DBGkmer& dbgkmer) {
        std::copy(std::begin(dbgkmer.fw), std::end(dbgkmer.fw), std::begin(fw));
        std::copy(std::begin(dbgkmer.bw), std::end(dbgkmer.bw), std::begin(bw));
        cov = dbgkmer.cov;
    }
    
    uint8_t fwCount() { // count forward edges
        
        uint32_t fwEdges = 0;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (fw[i] != 0)
                ++fwEdges;
        }
        return fwEdges;
    }
    
    uint8_t bwCount() { // count backward edges
        
        uint32_t bwEdges = 0;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (bw[i] != 0)
                ++bwEdges;
        }
        return bwEdges;
    }
    
    std::vector<uint8_t> fwEdgeIndexes() { // count forward edges
        
        std::vector<uint8_t> indexes;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (fw[i] != 0)
                indexes.push_back(i);
        }
        return indexes;
    }
    
    std::vector<uint8_t> bwEdgeIndexes() { // count forward edges
        
        std::vector<uint8_t> indexes;
        
        for (uint8_t i = 0; i<4; ++i) {
            if (bw[i] != 0)
                indexes.push_back(i);
        }
        return indexes;
    }
    
};

struct DBGkmer32color : public DBGkmer32 {
    uint8_t color = 0;
    
    DBGkmer32color() {}
    
    DBGkmer32color(const DBGkmer32& dbgkmer32) {
        std::copy(std::begin(dbgkmer32.fw), std::end(dbgkmer32.fw), std::begin(fw));
        std::copy(std::begin(dbgkmer32.bw), std::end(dbgkmer32.bw), std::begin(bw));
        cov = dbgkmer32.cov;
    }
};

template<typename T>
class PM : public phmap::parallel_flat_hash_map<uint64_t, T,
                                          std::hash<uint64_t>,
                                          std::equal_to<uint64_t>,
                                          std::allocator<std::pair<const uint64_t, T>>,
                                          8,
                                          phmap::NullMutex> {};

using ParallelMap = PM<DBGkmer>;
using ParallelMap32 = PM<DBGkmer32>;
using ParallelMap32color = PM<DBGkmer32color>;

class DBG : public Kmap<DBG, UserInputKreeq, uint64_t, DBGkmer, DBGkmer32> { // CRTP
    
    std::atomic<uint64_t> totMissingKmers{0}, totKcount{0}, totEdgeMissingKmers{0};
    UserInputKreeq userInput;
    
    InSequencesDBG *genome;
    
    ParallelMap32 *graphCache = new ParallelMap32;
    
    // subgraph objects
    ParallelMap32color *DBGsubgraph = new ParallelMap32color;
    std::vector<ParallelMap32color*> DBGTmpSubgraphs;
    InSequences GFAsubgraph;

    uint64_t totEdgeCount = 0;

public:
    
    DBG(UserInputKreeq& userInput) : Kmap{userInput}, userInput(userInput) {
        DBextension = "kreeq";
        
        if (userInput.kmerDepth == -1) {
            if (userInput.travAlgorithm == "best-first")
                this->userInput.kmerDepth = userInput.kmerLen; // unidirectional search
            else if (userInput.travAlgorithm == "traversal")
                this->userInput.kmerDepth = std::ceil((float)userInput.kmerLen/2); // kmer search is in both directions
        }
    }

    ~DBG(){
        delete DBGsubgraph;
    };
    
    bool hashSequences();
    
    bool processBuffers(uint16_t m);
    
    bool summary(uint16_t m);
    
    void DBstats();
    
    void loadGenome(InSequencesDBG *genome);
    
    void correctSequences();
    
    bool evaluateSegment(uint32_t i, std::array<uint16_t, 2> mapRange);
    
    void validateSequences();
    
    bool reloadMap32(uint16_t m);
    
    bool deleteMap(uint16_t m);
    
    bool mergeSubMaps(ParallelMap* map1, ParallelMap* map2, uint8_t subMapIndex, uint16_t m);
    
    void kunion();
    
    void report();
    
    void printTable(std::string ext);
    
    void printTableCompressed();
    
    void writeIndex(std::ofstream &ofs);
    
    void printTableCompressedBinary();
    
    void printVCF();
    
    void searchGraph();
    
    std::pair<DBGkmer*,bool> findDBGkmer(uint8_t *origin);
    
    int8_t scorePath(std::string path);
    
    int8_t checkNext(uint8_t *nextKmer, uint8_t *origin);
    
    std::deque<DBGpath> findPaths(uint8_t *origin, uint8_t *target, uint8_t depth, DBGpath currentPath, Log &threadLog);
    
    void printAltPaths(std::vector<std::vector<uint8_t>> altPaths, Log &threadLog);
    
    bool loadSegmentCoordinates(InSegment *inSegment, std::vector<uint64_t> &segmentCoordinates);
    
    BedCoordinates BEDPathsToSegments();
    
    bool detectAnomalies(InSegment *inSegment, std::vector<uint64_t> &anomalies);
    
    bool DBGtoVariants(InSegment *inSegment);
    
    std::pair<bool,std::deque<DBGpath>> searchVariants(std::pair<const uint64_t,DBGkmer32> source, bool isSourceFw, uint64_t ref, std::array<uint16_t, 2> mapRange, const std::deque<uint64_t> &targetsQueue, const phmap::parallel_flat_hash_map<uint64_t,bool> &targetsMap, ParallelMap32* localGraphCache);
    
    bool variantsToGFA(InSegment *inSegment, Log &threadLog);
    
    void collapseNodes();
    
    void DBGgraphToGFA();
    
    // subgraph functions
    
    bool DBGsubgraphFromSegment(InSegment *inSegment, std::array<uint16_t, 2> mapRange);
    
    void subgraph();
    
    void traversal(); // a basic traversal, all nodes at the same time
    
    void bestFirst(); // best-first search traversal
    
    void removeMissingEdges();
    
    void summary(ParallelMap32color& DBGsubgraph);
    
    ParallelMap32color traversalPass(ParallelMap32color* candidates, std::array<uint16_t, 2> mapRange);
    
    void mergeSubgraphs();
    
    std::pair<bool,ParallelMap32color> dijkstra(std::pair<uint64_t,DBGkmer32color>, std::array<uint16_t, 2> mapRange);
    
    void buildNextKmer(uint8_t* nextKmer, uint64_t hash, uint8_t nextBase, bool fw);
    
    void nextKmerFromString(uint8_t *nextKmer, std::string *sequence, uint64_t start, uint8_t nextBase);
    
    template<typename MAPTYPE>
    bool mergeSubMaps(MAPTYPE* map1, MAPTYPE* map2, uint8_t subMapIndex);
    
    template<typename MAPTYPE>
    bool unionSum(MAPTYPE* map1, MAPTYPE* map2);
    
};

#endif /* KREEQ_H */
