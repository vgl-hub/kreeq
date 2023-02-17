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

void Input::load(UserInputKreeq userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(bool mode, InSequences& inSequences) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    if (mode == 0) {
        
        DBG knav(userInput.kmerLen);
        
        lg.verbose("Loading input reads");
        
        unsigned int numFiles = userInput.iReadFileArg.size();
        
        for (unsigned int i = 0; i < numFiles; i++)
            loadSequences(userInput, &knav, 'r', &i);
        
        lg.verbose("Reads loaded");
        
        knav.build();
        
        knav.validateSequences(inSequences);
        
    }else{
        
    }
    
}

void Input::read(InSequences& inSequences) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    stream = streamObj.openStream(userInput, 'f');
    
    readGFA(inSequences, userInput, stream);

    jobWait(threadPool);
    
    inSequences.updateStats();

}
