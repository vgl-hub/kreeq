#include <stdlib.h>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>

#include <parallel_hashmap/phmap.h>
#include "parallel_hashmap/phmap_dump.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "global.h"
#include "uid-generator.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"
#include "input.h"

#include "kmer.h"
#include "kreeq.h"

bool Kpos::validate(UserInputKreeq userInput) {
    
    lg.verbose("Navigating with " + std::to_string(mapCount) + " maps");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return joinBuff(m); });
    
    jobWait(threadPool);
    
    lg.verbose("Generate summary statistics");
    
    for(uint16_t m = 0; m<mapCount; ++m)
        threadPool.queueJob([=]{ return stats(map[m]); });
    
    jobWait(threadPool);
    
}

bool Kpos::joinBuff(uint16_t m) {
    
    
}

bool Kpos::stats(phmap::flat_hash_map<uint64_t, std::vector<pos>>& map) {
    
    
}
