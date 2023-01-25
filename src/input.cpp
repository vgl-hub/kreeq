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
#include "sak.h" // swiss army knife
#include "zlib.h"
#include "stream-obj.h"
#include "input-gfa.h"

#include "kreeq.h"
#include "kmer.h"
#include "input.h"

void Input::load(UserInputKreeq userInput) {
    
    this->userInput = userInput;
    
}

void Input::read(bool mode) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    if (mode == 0) {
        
        Kpos knav(userInput.kmerLen);
        
        lg.verbose("Kmer object generated");
        
        knav.convert(userInput);
        
        knav.validate(userInput);
        
        knav.report(userInput);
        
    }else{
        
    }
    
}
