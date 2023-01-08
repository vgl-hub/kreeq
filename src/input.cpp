#include <stdlib.h>
#include <unistd.h>
#include <string>
#include <thread>
#include <mutex>
#include <vector>
#include <queue>
#include <stack>
#include <math.h>

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

#include "ktree.h"
#include "input.h"

inline size_t Input::hash(unsigned short int *kmer) {
    size_t result = 0;
    for(unsigned short int c = 0; c<k; c++) {
        
        result += *kmer++ * pow(4,c);
    }
    return result;
}

bool Input::threadedInsert(std::vector<InSegment*>& segments, unsigned short int thread_idx) {
    
    unsigned short int mapCount = 16;

    unsigned long long int moduloMap = (unsigned long long int) pow(4,k) / mapCount;
    
    unsigned short int moduloThread = mapCount / (threadPool.totalThreads() - 1);
    
    for (InSegment* segment : segments) {
        
        if (segment->getSegmentLen()<k)
            continue;

        unsigned long long int len = segment->getSegmentLen()-k+1;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        unsigned short int* str = new unsigned short int[segment->getSegmentLen()];
        
        for (unsigned long long int i = 0; i < len+k-1; i++){
            
            str[i] = ctoi[*(first+i)];
            
        }
        
        for (unsigned long long int c = 0; c<len; ++c){
            
            unsigned long long int value = hash(str+c);
            
            unsigned long long int map = value / moduloMap;
            
            if (map / moduloThread == thread_idx) // if the key belongs to this thread
                ++kcount[map][value];
            
//            std::cout<<value<<" "<<map<<" "<<moduloMap<<" "<<moduloThread<<" "<<map / moduloThread <<"="<<thread_idx<<std::endl;
            
        }
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
        delete[] str;

    }
    
    return true;

}

void Input::load(UserInputKreeq userInput) {
    
    this->userInput = userInput;
    
}


void Input::read(InSequences& inSequences) {
    
    if (userInput.iSeqFileArg.empty()) {return;}
    
    threadPool.init(maxThreads); // initialize threadpool
    
    stream = streamObj.openStream(userInput, 'f');
    
    if (!userInput.iSeqFileArg.empty() || userInput.pipeType == 'f') {
        
        StreamObj streamObj;
        
        stream = streamObj.openStream(userInput, 'f');
        
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
                        
                        Sequence* sequence = new Sequence{seqHeader, seqComment, inSequence, NULL};
                        
                        sequence->seqPos = seqPos; // remember the order
                        
                        inSequences.appendSequence(sequence);
                        
                        seqPos++;
                        
                    }
                    
                    break;
                    
                }
                    
            }
            
        }
        
    }
    
    jobWait(threadPool);
    
    if(verbose_flag) {std::cerr<<"\n\n";};
    
    std::vector<Log> logs = inSequences.getLogs();
    
    //consolidate log
    for (auto it = logs.begin(); it != logs.end(); it++) {
        
        it->print();
        logs.erase(it--);
        if(verbose_flag) {std::cerr<<"\n";};
        
    }
    
    std::vector<InSegment*>* segments = inSequences.getInSegments();
    
    inSequences.updateStats();
    
    k = userInput.kmerLen;
    
    for(unsigned short int t = 0; t<threadPool.totalThreads(); t++)
        threadPool.queueJob([=]{ return threadedInsert(*segments, t); });
    
    jobWait(threadPool);
    
    unsigned long long int totKmersUnique = 0;
    
    for (unsigned short int m = 0; m < 16; m++) {
        
        std::cout<<"map: "<<m<<std::endl;
        print_map(kcount[m]);

        for (auto const& pair : kcount[m]) {

            if(pair.second > 0){
                totKmers += pair.second;
                ++totKmersUnique;
            }

        }

    }

    std::cout<<"Total kmers: "<<totKmers<<std::endl;
    std::cout<<"Unique kmers: "<<totKmersUnique<<std::endl;

//	print_map(kcount);
	
	threadPool.join();
    
}
