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
    
};

using parallelMap = phmap::parallel_flat_hash_map<uint64_t, DBGkmer,
                                          std::hash<uint64_t>,
                                          std::equal_to<uint64_t>,
                                          std::allocator<std::pair<const uint64_t, DBGkmer>>,
                                          8,
                                          phmap::NullMutex>;

using parallelMap32 = phmap::parallel_flat_hash_map<uint64_t, DBGkmer32,
                                          std::hash<uint64_t>,
                                          std::equal_to<uint64_t>,
                                          std::allocator<std::pair<const uint64_t, DBGkmer32>>,
                                          8,
                                          phmap::NullMutex>;

class DBG : public Kmap<DBG, UserInputKreeq, DBGkmer, DBGkmer32> { // CRTP
    
    std::atomic<uint64_t> totMissingKmers{0}, totKcount{0}, totEdgeMissingKmers{0}, buffers{0};
    std::atomic<bool> readingDone{false};
    std::vector<std::thread> threads;
    std::vector<std::future<bool>> futures;
    UserInputKreeq userInput;
    
    InSequencesDBG *genome;
    
    // subgraph objects
    parallelMap32 *DBGsubgraph = new parallelMap32;
    std::vector<parallelMap32*> DBGTmpSubgraphs;
    InSequences GFAsubgraph;
    
    std::queue<std::string*> readBatches;
    uint64_t totEdgeCount = 0;

public:
    
    DBG(UserInputKreeq& userInput) : Kmap{userInput}, userInput(userInput) {
        DBextension = "kreeq";
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
    
    bool mergeSubMaps(parallelMap* map1, parallelMap* map2, uint8_t subMapIndex, uint16_t m);
    
    bool mergeSubMaps(parallelMap32* map1, parallelMap32* map2, uint8_t subMapIndex);
    
    bool unionSum(parallelMap32* map1, parallelMap32* map2);
    
    void kunion();
    
    void report();
    
    void printTable(std::string ext);
    
    void printTableCompressed();
    
    void writeIndex(std::ofstream &ofs);
    
    void printTableCompressedBinary();
    
    void printVCF();
    
    bool searchGraph(std::array<uint16_t, 2> mapRange);
    
    std::pair<DBGkmer*,bool> findDBGkmer(uint8_t *origin);
    
    int8_t scorePath(std::string path);
    
    int8_t checkNext(uint8_t *nextKmer, uint8_t *origin);
    
    std::deque<DBGpath> findPaths(uint8_t *origin, uint8_t *target, uint8_t depth, DBGpath currentPath, Log &threadLog);
    
    void printAltPaths(std::vector<std::vector<uint8_t>> altPaths, Log &threadLog);
    
    bool loadSegmentCoordinates(InSegment *inSegment, std::vector<uint64_t> &segmentCoordinates);
    
    BedCoordinates BEDPathsToSegments();
    
    bool detectAnomalies(InSegment *inSegment, std::vector<uint64_t> &anomalies);
    
    bool DBGtoVariants(InSegment *inSegment);
    
    bool variantsToGFA(InSegment *inSegment, Log &threadLog);
    
    bool DBGsubgraphFromSegment(InSegment *inSegment, std::array<uint16_t, 2> mapRange);
    
    void subgraph();
    
    void mergeSubgraphs();
    
    void DBGgraphToGFA();
    
};

#endif /* KREEQ_H */
