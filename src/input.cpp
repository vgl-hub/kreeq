#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>

#include <iostream>
#include <fstream>
#include <sstream>

#include <parallel_hashmap/phmap.h>

#include "bed.h"
#include "struct.h"
#include "functions.h" // global functions
#include "log.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "threadpool.h"
#include "gfa.h"
#include "stream-obj.h"
#include "input-gfa.h"
#include "fastx.h"

#include "input.h"
#include "kmer.h"
#include "kreeq.h"

void Input::loadInput(UserInputKreeq userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(bool mode, InSequences& inSequences) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    if (mode == 0) { // sequence validation
        
        DBG knav(userInput.kmerLen); // navigational kmerdb
        
        if (userInput.iReadFileArg.size() > 0) {
            
            lg.verbose("Loading input reads");
            
            unsigned int numFiles = userInput.iReadFileArg.size(); // number of input files
            
            for (unsigned int i = 0; i < numFiles; i++) // load each input file in the kmerdb
                loadKmers(userInput, &knav, 'r', &i);
            
            lg.verbose("Reads loaded");
            
        }else{
            
            std::ifstream file;
            
            file.open(userInput.iSeqFileArg + "/.index"); // reads the kmer length from the index file
            std::string line;
            
            getline(file, line);
            file.close();
            
            knav.load(userInput); // loads kmers into the new kreeqdb
            
        }
        
        knav.finalize(); // populate the hash table
        
        knav.validateSequences(inSequences); // validate the input sequence
        
        knav.report(userInput); // output
        
    }else{
        
    }
    
}

void Input::loadSequences(InSequences& inSequences) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    //intermediates
    std::string h;
    char* c;
    
    // stream read variable definition
    std::string firstLine;
    unsigned int seqPos = 0; // to keep track of the original sequence order
    
    std::string newLine, seqHeader, seqComment, line, bedHeader;
    
    stream = streamObj.openStream(userInput, 'f'); // open file
    
    if (stream) {
        
        switch (stream->peek()) {
                
            case '>': {
                
                stream->get();
                
                while (getline(*stream, newLine)) {
                    
                    h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                    c = strtok(NULL,""); //read comment
                    
                    seqHeader = h;
                    
                    if (c != NULL) {
                        
                        seqComment = std::string(c);
                        
                    }
                    
                    std::string* inSequence = new std::string;
                    
                    getline(*stream, *inSequence, '>');
                    
                    lg.verbose("Individual fasta sequence read");
                    
                    Sequence* sequence = new Sequence {seqHeader, seqComment, inSequence};
                    
                    if (sequence != NULL) {
                        
                        sequence->seqPos = seqPos; // remember the order
                        
                        inSequences.appendSequence(sequence);
                        
                        seqPos++;
                        
                    }
                    
                }
                
                break;
            }
            case '@': {
                
                while (getline(*stream, newLine)) { // file input
                    
                    newLine.erase(0, 1);
                    
                    h = std::string(strtok(strdup(newLine.c_str())," ")); //process header line
                    c = strtok(NULL,""); //read comment
                    
                    seqHeader = h;
                    
                    if (c != NULL) {
                        
                        seqComment = std::string(c);
                        
                    }
                    
                    std::string* inSequence = new std::string;
                    getline(*stream, *inSequence);
                    
                    getline(*stream, newLine);
                    
                    std::string* inSequenceQuality = new std::string;
                    getline(*stream, *inSequenceQuality);
                    
                    Sequence* sequence = new Sequence {seqHeader, seqComment, inSequence, inSequenceQuality};
                    
                    if (sequence != NULL) {
                        
                        sequence->seqPos = seqPos; // remember the order
                    
                        inSequences.appendSequence(sequence);
                        
                        seqPos++;
                        
                    }
                    
                }
                
                break;
                
            }
            default: {
                
                readGFA(inSequences, userInput, stream);
                
            }
            
        }
        
        lg.verbose("End of file");
            
    }else{

        fprintf(stderr, "Stream not successful: %s", userInput.iSeqFileArg.c_str());
        exit(1);

    }

    jobWait(threadPool);
    
    inSequences.updateStats(); // compute summary statistics

}
