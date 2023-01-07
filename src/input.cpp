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

inline size_t Input::hash(unsigned short int *kmer)
{
    size_t result = 0;
    for(unsigned short int c = 0; c<k; c++) {
        
        result += *kmer+c * pow(4,c);
    }
    return result;
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
    
    unsigned int* kcount = new unsigned int[(unsigned long long int)pow(4,k)]();

	for (InSegment* segment : *segments) {
        
        if (segment->getSegmentLen()<k)
            continue;

        unsigned long long int len = segment->getSegmentLen()-k+1;
        
        unsigned char* first = (unsigned char*)segment->getInSequencePtr()->c_str();
        
        unsigned short int* str = new unsigned short int[segment->getSegmentLen()];
        
        for (unsigned long long int i = 0; i < len+k-1; i++){
            
            str[i] = ctoi[*(first+i)];
            
        }

        for (unsigned long long int c = 0; c<len; ++c){
            
            ++kcount[hash(str+c)];
            
        }
        
        lg.verbose("Processed segment: " + segment->getSeqHeader());
        
        delete[] str;

	}
    
    unsigned long long int totKmersUnique = 0;
    
    for (unsigned long long int c = 0; c<pow(4,k); ++c) {
        
        if(kcount[c] > 0){
            totKmers += kcount[c];
            ++totKmersUnique;
        }

    }

    std::cout<<"Total kmers: "<<totKmers<<std::endl;
    std::cout<<"Unique kmers: "<<totKmersUnique<<std::endl;
    
    delete[] kcount;

//	print_map(kcount);
	
//	Ktree ktree(inSequences, userInput.kmerLen);
	
//	std::vector<InSegment*>* segments = inSequences.getInSegments();
//
//	std::vector<Knode*> nodes;
//
//	const short int k = userInput.kmerLen;
//
//	unsigned long long i = 0;
//
//	for (InSegment* segment : *segments) {
//
//		lg.verbose("Processing segment: " + segment->getSeqHeader());
//
//		unsigned long long int len = segment->getSegmentLen()-k+1;
//
//		std::string* seq = segment->getInSequencePtr();
//
//		char* first = segment->first();
//
//		for (unsigned long long int c = 0; c<len; ++c) {
//
//			std::string kmer_for = seq->substr(c, k);
//
//			Knode* fw = new Knode(first+c);
//
//			Knode* rc = new Knode(first+c);
//
//			for (unsigned short int j = 0; j<k; ++j) {
//
//				++i;
//
//			}
//
//			delete fw, delete rc;
//
//		}
//
//	}
//
//	std::cout<<nodes.size()<<" "<<i<<std::endl;
	
	threadPool.join();
    
}
