#ifndef INPUT_H
#define INPUT_H

struct UserInputKreeq : UserInput {

    unsigned short int kmerLen = 21;

};

class Input {
    
    UserInputKreeq userInput;
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    StreamObj streamObj;
    
    std::shared_ptr<std::istream> stream;
    
    unsigned int totKmers = 0;
    
public:
    
    void load(UserInputKreeq userInput);
    
    void read(InSequences& inSequence);
    
    inline size_t hash(const char * string);
    
};

#endif /* INPUT_H */
