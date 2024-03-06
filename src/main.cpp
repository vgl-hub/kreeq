#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <string>
#include <getopt.h>
#include <fstream>

#include "global.h"
#include "log.h"
#include "uid-generator.h"

#include "parallel-hashmap/phmap.h"

#include "bed.h"
#include "struct.h"
#include "functions.h"
#include "gfa-lines.h"
#include "gfa.h"
#include "stream-obj.h"

#include "input.h"
#include "main.h"

std::string version = "0.0.1";

//global
std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run

int cmd_flag;
int verbose_flag;
short int tabular_flag;
int maxThreads;

std::mutex mtx;
ThreadPool<std::function<bool()>> threadPool;
Log lg;
std::vector<Log> logs;

void printHelp() {
    
    printf("kreeq [mode]\n-h for additional help.\n");
    printf("\nModes:\n");
    printf("validate\n");
    printf("union\n");
    exit(0);
    
}

int main(int argc, char **argv) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    bool isPipe = false; // to check if input is from pipe
    uint8_t mode = 0;
    
    std::string cmd;

    UserInputKreeq userInput; // initialize input object
    
    if (argc == 1) { // gfastats with no arguments
            
        printHelp();
        
    }

    if(strcmp(argv[1],"validate") == 0) {

        mode = 0;
        
        static struct option long_options[] = { // struct mapping long options
            {"coverage-cutoff", required_argument, 0, 'c'},
            {"database", required_argument, 0, 'd'},
            {"input-sequence", required_argument, 0, 'f'},
            {"kmer-length", required_argument, 0, 'k'},
            {"search-depth", required_argument, 0, 0},
            {"backtracking-span", required_argument, 0, 0},
            {"out-format", required_argument, 0, 'o'},
            {"input-reads", required_argument, 0, 'r'},
            {"tmp-prefix", required_argument, 0, 't'},
            {"max-memory", required_argument, 0, 'm'},
            {"threads", required_argument, 0, 'j'},
            
            {"verbose", no_argument, &verbose_flag, 1},
            {"cmd", no_argument, &cmd_flag, 1},
            {"version", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            
            
            {0, 0, 0, 0}
        };
        
        while (arguments) { // loop through argv
            
            int option_index = 0;
            
            c = getopt_long(argc, argv, "-:c:d:f:k:o:r:t:m:j:v:h",
                            long_options, &option_index);
            
            if (c == -1) { // exit the loop if run out of options
                break;
                
            }
            
            switch (c) {
                case ':': // handle options without arguments
                    switch (optopt) { // the command line option last matched
                        case 'b':
                            break;
                            
                        default:
                            fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                            return EXIT_FAILURE;
                    }
                    break;
                default: // handle positional arguments
                    
                case 0: // case for long options without short options
                    
                    if(strcmp(long_options[option_index].name,"search-depth") == 0)
                        userInput.depth = atoi(optarg);
                    
                    if(strcmp(long_options[option_index].name,"backtracking-span") == 0)
                        userInput.backtrackingSpan = atoi(optarg);
                    
                    break;

                case 'c': // input kreeq db
                    
                    if (!isNumber(optarg)) {
                        fprintf(stderr, "input '%s' to option -%c must be a number\n", optarg, optopt);
                        return EXIT_FAILURE;
                    }
                    
                    userInput.covCutOff = atoi(optarg);
                    
                    break;

                case 'd': // input kreeq db
                    
                    if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                        
                        userInput.pipeType = 'f'; // pipe input is a sequence
                        
                    }else{ // input is a regular file
                        
                        ifFileExists(optarg);
                        userInput.inDBG.push_back(optarg);
                        
                    }
                    
                    break;
                    
                case 'f': // input sequence
                    
                    if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                        
                        userInput.pipeType = 'f'; // pipe input is a sequence
                        
                    }else{ // input is a regular file
                        
                        ifFileExists(optarg);
                        userInput.inSequence = optarg;
                        
                    }
                    
                    break;
                    
                case 'k': // kmer length
                    
                    if (!isNumber(optarg)) {
                        fprintf(stderr, "input '%s' to option -%c must be a number\n", optarg, optopt);
                        return EXIT_FAILURE;
                    }
                    
                    userInput.kmerLen = atoi(optarg);
                    
                    break;
                    
                case 'j': // max threads
                    maxThreads = atoi(optarg);
                    break;

                case 'o': // handle output (file or stdout)
                    userInput.outFile = optarg;
                    break;

                case 'r': // input reads
                    
                    if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                    
                        userInput.pipeType = 'r'; // pipe input is a sequence
                    
                    }else{ // input is a regular file
                        
                        optind--;
                        for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
                            
                            ifFileExists(argv[optind]);
                            userInput.inReads.push_back(argv[optind]);
                            
                        }
                    }
                        
                    break;
                    
                case 't': // prefix for temporary files
                    userInput.prefix = optarg;
                    break;
                    
                case 'm': // prefix for temporary files
                    userInput.maxMem = atof(optarg);
                    break;
                    
                case 'v': // software version
                    printf("kreeq v%s\n", version.c_str());
                    printf("Giulio Formenti giulio.formenti@gmail.com\n");
                    exit(0);
                    
                case 'h': // help
                    printf("kreeq [command]\n");
                    printf("\nOptions:\n");
                    printf("\t-d --database kreeq database to load.\n");
                    printf("\t-f --input-sequence sequence input file (fasta,gfa1/2).\n");
                    printf("\t-r --input-reads read input files (fastq).\n");
                    printf("\t-k --kmer-length length of kmers.\n");
                    printf("\t-o --out-format supported extensions:\n");
                    printf("\t\t .kreeq dumps hashmaps to file for reuse.\n");
                    printf("\t-t --tmp-prefix prefix to temporary directory.\n");
                    printf("\t-m --max-memory use at most this amount of memory (in Gb, default: 0.5 of max).\n");
                    printf("\t-j --threads <n> numbers of threads (default: max).\n");
                    printf("\t-v --version software version.\n");
                    printf("\t--cmd print $0 to stdout.\n");
                    exit(0);
            }
            
            if    (argc == 2 || // handle various cases in which the output should include default outputs
                   (argc == 3 && pos_op == 2) ||
                   (argc == 4 && pos_op == 3)) {
                
            }
            
        }
        
    }else if(strcmp(argv[1],"union") == 0){
        
        mode = 1;
        
        static struct option long_options[] = { // struct mapping long options
            {"databases", required_argument, 0, 'd'},
            {"out-format", required_argument, 0, 'o'},
            
            {"threads", required_argument, 0, 'j'},
            {"verbose", no_argument, &verbose_flag, 1},
            {"cmd", no_argument, &cmd_flag, 1},
            {"help", no_argument, 0, 'h'},
            
            {0, 0, 0, 0}
        };
        
        while (true) { // loop through argv
            
            int option_index = 1;
            
            c = getopt_long(argc, argv, "-:d:j:o:h",
                            long_options, &option_index);
            
            if (c == -1) { // exit the loop if run out of options
                break;
                
            }
            
            switch (c) {
                case ':': // handle options without arguments
                    switch (optopt) { // the command line option last matched
                        case 'b':
                            break;
                            
                        default:
                            fprintf(stderr, "option -%c is missing a required argument\n", optopt);
                            return EXIT_FAILURE;
                    }
                    break;
                default: // handle positional arguments
                    
                case 0: // case for long options without short options
                    
                    //                if (strcmp(long_options[option_index].name,"line-length") == 0)
                    //                  splitLength = atoi(optarg);
                    
                    break;
                    
                case 'd': // input sequence
                    
                    if (isPipe && userInput.pipeType == 'n') { // check whether input is from pipe and that pipe input was not already set
                        
                        userInput.pipeType = 'f'; // pipe input is a sequence
                        
                    }else{ // input is a regular file
                        
                        optind--;
                        for( ;optind < argc && *argv[optind] != '-' && !isInt(argv[optind]); optind++){
                            
                            ifFileExists(argv[optind]);
                            userInput.inDBG.push_back(argv[optind]);
                            
                        }
                        
                    }
                    
                    break;
                    
                case 'j': // max threads
                    maxThreads = atoi(optarg);
                    break;
                    
                case 'o': // handle output (file or stdout)
                    userInput.outFile = optarg;
                    break;
                    
                case 'h': // help
                    printf("kreeq union [options]\n");
                    printf("\nOptions:\n");
                    printf("\t-d --databases DBG databases to merge.\n");
                    printf("\t-j --threads <n> numbers of threads (default: max).\n");
                    printf("\t-o --out-format generates various kinds of outputs (currently supported: .kreeq).\n");
                    printf("\t--cmd print $0 to stdout.\n");
                    exit(0);
            }
            
            if    (argc == 2 || // handle various cases in which the output should include default outputs
                   (argc == 3 && pos_op == 2) ||
                   (argc == 4 && pos_op == 3)) {
                
            }
            
        }
        
        if (userInput.inDBG.size() < 2) {
            fprintf(stderr, "At least two databases required (-d).\n");
            return EXIT_FAILURE;
        }
        
    }else{
        
        fprintf(stderr, "Unrecognized mode: %s\n", argv[1]);
        printHelp();
        
    }
    
    if (cmd_flag) { // print command line
        for (unsigned short int arg_counter = 0; arg_counter < argc; arg_counter++) {
            printf("%s ", argv[arg_counter]);
        }
        printf("\n");
        
    }
    
    Input in;
    
    threadPool.init(maxThreads); // initialize threadpool
    maxMem = (userInput.maxMem == 0 ? get_mem_total(3) * 0.5 : userInput.maxMem); // set memory limit
    
    in.loadInput(userInput); // load user input
    lg.verbose("User input loaded");
    
    in.read(mode); // read input reads and validate

    threadPool.join(); // join threads
    
    exit(EXIT_SUCCESS);
	
}
