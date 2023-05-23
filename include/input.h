#ifndef INPUT_H
#define INPUT_H

struct UserInputKreeq : UserInput {

    uint8_t kmerLen = 21, covCutOff = 0;
    std::string outFile = "", iDBGFileArg = "";

};

class Input {
    
    UserInputKreeq userInput;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
        
public:
    
    void loadInput(UserInputKreeq userInput);
    
    void loadSequences(InSequences& inSequences);
    
    void read(bool mode, InSequences& inSequences);
    
};

#endif /* INPUT_H */
