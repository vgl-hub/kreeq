//
//  compressor.cpp
//  kreeq-dev
//
//  Created by Giulio Formenti on 2/11/24.
//

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <vector>
#include <regex>
#include <array>

int main(int argc, char *argv[])
{
    try
    {
        if (argc == 1) {
            std::cerr << "Please provide input file\n";
            return EXIT_FAILURE;
        }
            
        std::string fl = argv[1];
        std::ifstream ifs(fl, std::ios::in | std::ios::binary);

        uint8_t k;
        uint16_t size;
        uint64_t len;
        std::array<uint8_t, 3> values;
        char entrySep = ',', colSep = ',';
        uint64_t absPos = 0;
        
        ifs.read(reinterpret_cast<char *>(&k), sizeof(uint8_t));
        if (argv[2] == NULL)
            std::cout<<std::to_string(k)<<"\n";

        while(ifs && !(ifs.peek() == EOF)) {
            
            ifs.read(reinterpret_cast<char *>(&size), sizeof(uint16_t));
            std::string header;
            header.resize(size);
            ifs.read(reinterpret_cast<char *>(&header[0]), sizeof(header[0])*size);
            ifs.read(reinterpret_cast<char *>(&len), sizeof(uint64_t));

            std::regex start("start\\=(\\d+)");
            std::smatch match;
            std::regex_search(header, match, start);
            absPos = std::stoi(match[1]);
            
            if(argv[2] == NULL) {
                
                std::cout<<header<<"\n";
                    
                for (uint64_t i = 0; i < len; ++i) {
                    ifs.read(reinterpret_cast<char *>(&values), sizeof(uint8_t)*3);
                    for (auto iter = values.begin(); iter != values.end(); iter++) {
                        std::cout<<std::to_string(*iter);
                        if (std::next(iter) != values.end())
                            std::cout<<",";
                    }
                    std::cout<<"\n";
                }
                    
            } else if (argv[2] == std::string("--expand")){
                
                std::regex chrom("chrom\\=(\\.*) ");
                std::regex_search(header, match, chrom);
                header = match[1];
                
                for (uint64_t i = 0; i < len; ++i) {
                    ifs.read(reinterpret_cast<char *>(&values), sizeof(uint8_t)*3);
                    
                    std::vector<uint8_t> kmerCov(k-1,0);
                    std::vector<uint8_t> edgeCovFw(k-1,0);
                    std::vector<uint8_t> edgeCovBw(k-1,0);
                    
                    std::cout<<header<<colSep<<absPos<<colSep;
                    
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
    catch(std::exception const& e)
    {
        std::cerr << e.what() << '\n';
        return EXIT_FAILURE;
    }

    return EXIT_SUCCESS;
}
