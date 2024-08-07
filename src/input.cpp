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

#include <parallel-hashmap/phmap.h>

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

InSequencesDBG::~InSequencesDBG() {
    for (DBGbase *p : dbgbases)
        delete p;
}

void InSequencesDBG::generateValidationVector() {

    for (InSegment* segment : inSegments) {
        alloc += segment->getSegmentLen() * sizeof(DBGbase);
        DBGbase *dbgbase = new DBGbase[segment->getSegmentLen()];
        dbgbases.push_back(dbgbase);
    }
}

std::vector<DBGbase*>* InSequencesDBG::getInSegmentsDBG() {
    return &dbgbases;
}

void Input::loadInput(UserInputKreeq userInput) {
    this->userInput = userInput;
}


void Input::loadGraph() {
    
    if (userInput.kmerDB.size() == 1){
        userInput.prefix = userInput.kmerDB[0]; // access database
        std::ifstream file;
        file.open(userInput.prefix + "/.index"); // update kmer length
        std::string line;
        getline(file, line);
        file.close();
        lg.verbose("Overriding default kmer length (" + std::to_string(userInput.kmerLen) + ") with DB kmer length (" + line + ").");
        userInput.kmerLen = stoi(line);
    }else if (userInput.kmerDB.size() > 1) {
        fprintf(stderr, "More than one DBG database provided. Merge them first. Exiting.\n");
        exit(EXIT_FAILURE);
    }else{
        fprintf(stderr, "Cannot load DBG input. Exiting.\n");
        exit(EXIT_FAILURE);
    }
}

void Input::read() {
    
    if (userInput.outFile.find(".kreeq") != std::string::npos)
        userInput.prefix = userInput.outFile;
    
    if (userInput.prefix != ".")
        make_dir(userInput.prefix.c_str());
    
    switch (userInput.mode) {
            
        case 0: { // sequence validation
            
            DBG knav(userInput); // navigational kmerdb

            if (userInput.inReads.size() > 0) {
                
                lg.verbose("Loading input reads.");
                unsigned int numFiles = userInput.inReads.size(); // number of input files
                
                for (unsigned int i = 0; i < numFiles; ++i) // load each input file in the kmerdb
                    loadKmers(userInput, knav, 'r', i);
                
                lg.verbose("Reads loaded.");
                knav.finalize();
            }else{
                loadGraph();
            }
            knav.loadHighCopyKmers(); // reload high copy kmers for computation steps
            
            InSequencesDBG genome; // initialize sequence collection object
            if (!userInput.inSequence.empty()) {
                
                if (!userInput.inSequence.empty()) {
                    lg.verbose("Loading input sequences");
                    loadGenome(genome); // read input genome
                    lg.verbose("Sequences loaded");
                }
                knav.loadGenome(&genome);
            }
            knav.report(); // output
            knav.cleanup(); // delete tmp files
            break;
        }
        case 1: { // union of multiple kmerdbs
        
            std::ifstream file;
            lg.verbose("Merging input databases.");
            unsigned int numFiles = userInput.kmerDB.size(); // number of input kmerdbs
            short unsigned int k = 0;
            
            for (unsigned int i = 0; i < numFiles; i++) {  // reads the kmer length from the index files checking consistency between kmerdbs
                
                file.open(userInput.kmerDB[i] + "/.index");
                std::string line;
                
                getline(file, line);
                file.close();
                
                if (k == 0)
                    k = stoi(line);
                
                if (k != stoi(line)) {
                    fprintf(stderr, "Cannot merge databases with different kmer length.\n");
                    exit(1);
                }
            }
            if (k == 0 || k > 32) {
                fprintf(stderr, "Invalid kmer length.\n");
                exit(1);
            }
            DBG knav(userInput); // a new empty kmerdb with the specified kmer length
            lg.verbose("DBG object generated. Merging.");
            knav.kunion(); // union set
            knav.report(); // output
            knav.cleanup(); // delete tmp files
            break;
        }
        case 2: { // subgraph
            
            loadGraph();
            DBG knav(userInput); // navigational kmerdb
            knav.loadHighCopyKmers(); // reload high copy kmers for computation steps

            InSequencesDBG genome; // initialize sequence collection object
            if (!userInput.inSequence.empty()) {
                
                if (!userInput.inSequence.empty()) {
                    lg.verbose("Loading input sequences");
                    loadGenome(genome); // read input genome
                    lg.verbose("Sequences loaded");
                }
                knav.loadGenome(&genome);
            }
            knav.subgraph();
            
            lg.verbose("Searching graph");
            knav.searchGraph();
            lg.verbose("Remove missing edges");
            knav.removeMissingEdges();
            lg.verbose("Generating GFA");
            knav.DBGgraphToGFA();
            
            knav.report(); // output
            knav.cleanup(); // delete tmp files
            break;
        }
        default:
            fprintf(stderr, "Invalid mode.\n");
            exit(1);
    }
}

void Input::loadGenome(InSequencesDBG& inSequences) {
    
    if (userInput.inSequence.empty()) {return;}
    
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
                    
                    if (c != NULL) 
                        seqComment = std::string(c);
                    
                    std::string* inSequence = new std::string;
                    
                    getline(*stream, *inSequence, '>');
                    alloc += inSequence->size() * sizeof(char);
                    
                    lg.verbose("Individual fasta sequence read");
                    
                    Sequence* sequence = new Sequence {seqHeader, seqComment, inSequence};
                    
                    if (sequence != NULL) {
                        
                        sequence->seqPos = seqPos; // remember the order
                        
                        inSequences.appendSequence(sequence, -1);
                        
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
                    
                        inSequences.appendSequence(sequence, -1);
                        
                        seqPos++;
                        
                    }
                    
                }
                
                break;
                
            }
            default: {
                
                readGFA(inSequences, userInput, stream);
                
            }
            
        }
        
        lg.verbose("End of file.");
            
    }else{

        fprintf(stderr, "Stream not successful: %s.", userInput.inSequence.c_str());
        exit(1);

    }

    jobWait(threadPool);
    
    inSequences.updateStats(); // compute summary statistics

}
