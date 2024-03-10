//
//  compressor.cpp
//  kreeq-dev
//
//  Created by Giulio Formenti on 2/11/24.
//

#include <getopt.h>
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <array>
#include <array>
#include <sstream>

#include "parallel-hashmap/phmap.h"

#include "global.h"
#include "log.h"
#include "bed.h"
#include "struct.h"
#include "functions.h"

std::string version = "0.0.1";

std::chrono::high_resolution_clock::time_point start = std::chrono::high_resolution_clock::now(); // immediately start the clock when the program is run
int cmd_flag;
int verbose_flag;
int tabular_flag;
int maxThreads;

struct BkwigIndexComponent {
    
    ComponentType componentType;
    uint64_t bytePos, absPos, len;
    uint8_t step;
    
};

struct BkwigIndex {
    
    phmap::parallel_flat_hash_map<std::string,std::vector<BkwigIndexComponent>> paths;
    uint64_t indexByteSize = 0;
    
    void printIndex() {
        
        for (auto pair : paths) {
            std::cout<<pair.first<<std::endl;
            for (auto comp : pair.second)
                std::cout<<+comp.componentType<<"\t"<<+comp.absPos<<"\t"<<+comp.len<<"\t"<<+comp.step<<std::endl;
        }
    }
};

struct UserInputDecompressor : UserInput {

    std::string inputFile, coordinateFile, outFile;
    uint64_t maxMem = 0;
    int expand = 0;
    uint32_t span = 0;
    BkwigIndex bkwigIndex;
};


void printHelp() {
    
    printf("decompressor [mode]\n-h for additional help.\n");
    printf("\nModes:\n");
    printf("inflate\n");
    exit(0);
    
}

void readIndex(std::ifstream &ifs, BkwigIndex &bkwigIndex) { // reads: nPaths, and for each path reads size and path header and nComponents, and for each segment 1) type 2) absPos 3) segment length
    
    uint64_t bytePos = 0;
    
    uint32_t nPaths;
    ifs.read(reinterpret_cast<char *>(&nPaths), sizeof(uint32_t));
    bkwigIndex.indexByteSize += sizeof(uint32_t);
    
    for (uint32_t i = 0; i < nPaths; ++i) {
        
        uint16_t headerSize;
        ifs.read(reinterpret_cast<char *>(&headerSize), sizeof(uint16_t));
        bkwigIndex.indexByteSize += sizeof(uint16_t);
        std::string pHeader;
        pHeader.resize(headerSize);
        ifs.read(reinterpret_cast<char *>(&pHeader[0]), sizeof(pHeader[0])*headerSize);
        bkwigIndex.indexByteSize += sizeof(pHeader[0])*headerSize;
        uint32_t nComponents;
        ifs.read(reinterpret_cast<char *>(&nComponents), sizeof(uint32_t));
        bkwigIndex.indexByteSize += sizeof(uint32_t);
        
        std::vector<BkwigIndexComponent> componentsVec;
        
        for (uint32_t e = 0; e < nComponents; ++e) {
            
            bytePos += headerSize * sizeof(char) + sizeof(uint64_t);
            
            ComponentType componentType;
            ifs.read(reinterpret_cast<char *>(&componentType), sizeof(ComponentType));
            bkwigIndex.indexByteSize += sizeof(ComponentType);
            bytePos += sizeof(ComponentType);
            uint64_t absPos;
            ifs.read(reinterpret_cast<char *>(&absPos), sizeof(uint64_t));
            bkwigIndex.indexByteSize += sizeof(uint64_t);
            bytePos += sizeof(uint64_t)*absPos;
            uint64_t len;
            ifs.read(reinterpret_cast<char *>(&len), sizeof(uint64_t));
            bkwigIndex.indexByteSize += sizeof(uint64_t);
            uint8_t step;
            ifs.read(reinterpret_cast<char *>(&step), sizeof(uint8_t));
            if (componentType == SEGMENT)
                bytePos += sizeof(uint64_t)*len;
            
            componentsVec.push_back({componentType, bytePos, absPos, len, step});
        }
        bkwigIndex.paths[pHeader] = componentsVec;
    }
}

int main(int argc, char *argv[]) {
    
    short int c; // optarg
    short unsigned int pos_op = 1; // optional arguments
    
    bool arguments = true;
    uint8_t mode = 0;
    
    std::string cmd;
    std::string inputFile, coordinateFile;
    
    UserInputDecompressor userInput; // initialize input object
    
    if (argc == 1) // decompressor with no arguments
        printHelp();

    if(strcmp(argv[1],"inflate") == 0) {

        mode = 0;
        
        static struct option long_options[] = { // struct mapping long options
            {"input-file", required_argument, 0, 'i'},
            {"coordinate-file", required_argument, 0, 'c'},
            {"span", required_argument, 0, 's'},
            {"out-format", required_argument, 0, 'o'},
            {"max-memory", required_argument, 0, 'm'},
            {"threads", required_argument, 0, 'j'},
            {"expand", no_argument, &userInput.expand, 1},
            
            {"verbose", no_argument, &verbose_flag, 1},
            {"cmd", no_argument, &cmd_flag, 1},
            {"version", no_argument, 0, 'v'},
            {"help", no_argument, 0, 'h'},
            
            
            {0, 0, 0, 0}
        };
        
        while (arguments) { // loop through argv
            
            int option_index = 0;
            
            c = getopt_long(argc, argv, "-:i:c:s:o:m:j:v:h",
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
                    
                    break;
                    
                case 'i': // input file
                    ifFileExists(optarg);
                    userInput.inputFile = optarg;
                    break;
                    
                case 'c': // input coordinate file
                    ifFileExists(optarg);
                    userInput.coordinateFile = optarg;
                    break;
                    
                case 's': // max threads
                    userInput.span = atoi(optarg);
                    break;
                    
                case 'j': // max threads
                    maxThreads = atoi(optarg);
                    break;
                    
                case 'o': // handle output (file or stdout)
                    userInput.outFile = optarg;
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
                    printf("\t-i --input-file kreeq input.\n");
                    printf("\t-c --coordinate-file sequence coordinates to extract.\n");
                    printf("\t-s --span <int> print context before and after coordinate positions.\n");
                    printf("\t-o --out-format supported extensions:\n");
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
    
    switch (mode) {
            
        case 0: // inflate
            
        {
            
            BedCoordinates bedIncludeList;
            
            if (!userInput.inBedExclude.empty()) {
                
                std::string line, bedHeader;
                uint64_t begin, end;
                std::ifstream ifs(userInput.coordinateFile, std::ifstream::in);
                
                while (getline(ifs, line)) {
                    
                    std::istringstream iss(line);
                    iss >> bedHeader >> begin >> end;
                    
                    bedIncludeList.pushCoordinates(bedHeader, begin, end);
                    
                }
                
            }
            
            std::ifstream ifs(userInput.inputFile, std::ios::in | std::ios::binary);
            
            uint8_t k;
            uint16_t size;
            std::array<uint8_t, 3> values;
            char entrySep = ',', colSep = ',';
            uint64_t absPos = 0;
            
            ifs.read(reinterpret_cast<char *>(&k), sizeof(uint8_t)); // read k length
            if (!userInput.expand)
                std::cout<<std::to_string(k)<<"\n"; // print k length
            
            readIndex(ifs, userInput.bkwigIndex);
//            userInput.bkwigIndex.printIndex();
            
            while(ifs && !(ifs.peek() == EOF)) { // read the entire file
                
                ifs.read(reinterpret_cast<char *>(&size), sizeof(uint16_t)); // read header length
                std::string pHeader;
                pHeader.resize(size); // resize string to fit the header
                ifs.read(reinterpret_cast<char *>(&pHeader[0]), sizeof(pHeader[0])*size); // read header into string
                std::vector<BkwigIndexComponent> components = userInput.bkwigIndex.paths[pHeader];
                
                for(BkwigIndexComponent component : components) {
                    
                    uint64_t len = component.len;
                    
                    if(!userInput.expand) {
                        
                        std::cout<<"fixedStep chrom="<<pHeader<<" start="<<std::to_string(component.absPos)<<" step="<<std::to_string(component.step)<<std::endl;
                        
                        uint8_t *values = new uint8_t[len*3];
                        ifs.read(reinterpret_cast<char *>(values), sizeof(uint8_t)*len*3);
                        std::ostringstream os;
                        uint8_t comma = 0;
                        
                        for (uint64_t i = 0; i < len; ++i) { // loop through all position for this record
                            
                            os<<std::to_string(values[i]);
                            if (comma < 2) {
                                os<<",";
                                ++comma;
                            }else{
                                os<<"\n";
                                comma = 0;
                            }
                        }
                        delete[] values;
                        auto str = os.str();
                        std::cout.write(str.c_str(), static_cast<std::streamsize>(str.size()));
                        
                    }else{
                        
                        std::vector<uint8_t> kmerCov(k-1,0);
                        std::vector<uint8_t> edgeCovFw(k-1,0);
                        std::vector<uint8_t> edgeCovBw(k-1,0);
                        
                        absPos = component.absPos;
                        
                        for (uint64_t i = 0; i < len; ++i) { // loop through all position for this record
                            ifs.read(reinterpret_cast<char *>(&values), sizeof(uint8_t)*3);
                            
                            std::cout<<pHeader<<colSep<<absPos<<colSep;
                            
                            kmerCov.push_back(values[0]);
                            
                            for(uint8_t c = 0; c<k; ++c){
                                std::cout<<std::to_string(kmerCov[c]);
                                if (c < k - 1)
                                    std::cout<<entrySep;
                            }
                            
                            std::cout<<colSep;
                            
                            edgeCovFw.push_back(values[1]);
                            
                            for(uint8_t c = 0; c<k; ++c){
                                std::cout<<std::to_string(edgeCovFw[c]);
                                if (c < k - 1)
                                    std::cout<<entrySep;
                            }
                            
                            std::cout<<colSep;
                            
                            edgeCovBw.push_back(values[2]);
                            
                            for(uint8_t c = 0; c<k; ++c){
                                std::cout<<std::to_string(edgeCovBw[c]);
                                if (c < k - 1)
                                    std::cout<<entrySep;
                            }
                            
                            std::cout<<"\n";
                            
                            kmerCov.erase(kmerCov.begin());
                            edgeCovFw.erase(edgeCovFw.begin());
                            edgeCovBw.erase(edgeCovBw.begin());
                            
                            ++absPos;
                        }
                    }
                }
            }
        }
    }
    return EXIT_SUCCESS;
}
