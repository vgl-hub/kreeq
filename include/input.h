#ifndef INPUT_H
#define INPUT_H

struct UserInputKreeq : UserInput {

    uint8_t kmerLen = 21;
    std::string outFile = "";

};

class Input {
    
    UserInputKreeq userInput;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
        
public:
    
    void load(UserInputKreeq userInput);
    
    void read(bool mode, InSequences& inSequences);
    
    void read(InSequences& inSequences);
    
};

#endif /* INPUT_H */
