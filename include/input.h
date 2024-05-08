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

    uint8_t kmerLen = 21, covCutOff = 0, depth = 3, backtrackingSpan = 5;
    uint64_t maxMem = 0;
    std::string prefix = ".", outFile = "";
    std::vector<std::string> inDBG;

};

class Input {
    
    UserInputKreeq userInput;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
        
public:
    
    void loadInput(UserInputKreeq userInput);
    
    void loadGenome(InSequencesDBG& inSequences);
    
    void read(uint8_t mode);
    
};

#endif /* INPUT_H */
