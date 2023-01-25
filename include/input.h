#ifndef INPUT_H
#define INPUT_H

struct UserInputKreeq : UserInput {

    unsigned short int kmerLen = 21;
    std::string outFile = "";

};

class Input {
    
    UserInputKreeq userInput;
    
    //intermediates
    std::string h;
    
    // stream read variable definition
    std::string firstLine;
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
    
public:
    
    void load(UserInputKreeq userInput);
    
    void read(bool mode);
    
};

#endif /* INPUT_H */
