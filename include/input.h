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

    uint32_t covCutOff = 0;
    int16_t kmerDepth = -1;
    uint8_t maxSpan = 5;
    // bool
    int noCollapse = 0, noReference = 0;
    std::string travAlgorithm = "best-first";

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
