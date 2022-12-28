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

#include "ktree.h"
#include "input.h"

void Input::load(UserInput userInput) {
    
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
	
	std::cout<<(*segments)[0]->getInSequence()<<std::endl;
	
	threadPool.join();
    
}
