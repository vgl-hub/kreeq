#ifndef INPUT_H
#define INPUT_H

struct DBGbase {
    
    uint32_t fw = 0, bw = 0, cov = 0;
    bool isFw = false;
    
};

class InSequencesDBG : public InSequences {
    
    std::vector<DBGbase*> dbgbases;
    
public:
    
    ~InSequencesDBG();
    
    void generateValidationVector();
    
    std::vector<DBGbase*>* getInSegmentsDBG();
    
};

struct UserInputKreeq : UserInput {

    uint32_t covCutOff = 0, kmerDepth = std::ceil((float)kmerLen/2); // kmer search is in both directions
    uint8_t depth = 3, backtrackingSpan = 5;
    uint64_t maxMem = 0;
    // bool
    int noCollapse = 0, noReference = 0;

};

class Input {
    
    UserInputKreeq userInput;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
        
public:
    
    void loadInput(UserInputKreeq userInput);
    
    void loadGenome(InSequencesDBG& inSequences);
    
    void loadGraph();
    
    void read();
    
};

#endif /* INPUT_H */
