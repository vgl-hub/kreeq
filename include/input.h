#ifndef INPUT_H
#define INPUT_H

struct UserInputKreeq : UserInput {

    uint8_t kmerLen = 21;
    std::string outFile = "";

};

class Input {
    
    UserInputKreeq userInput;
        
public:
    
    void load(UserInputKreeq userInput);
    
    void read(bool mode);
    
};

#endif /* INPUT_H */
