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
#include <array>
#include <sstream>

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
        
        ifs.read(reinterpret_cast<char *>(&k), sizeof(uint8_t)); // read k length
        if (argv[2] == NULL)
            std::cout<<std::to_string(k)<<"\n"; // print k length

        while(ifs && !(ifs.peek() == EOF)) { // read the entire file
            
            ifs.read(reinterpret_cast<char *>(&size), sizeof(uint16_t)); // read header length
            std::string header;
            header.resize(size); // resize string to fit the header
            ifs.read(reinterpret_cast<char *>(&header[0]), sizeof(header[0])*size); // read header into string
            ifs.read(reinterpret_cast<char *>(&len), sizeof(uint64_t)); // read number of bases/rows

            std::regex start("start\\=(\\d+)"); // match wig start position
            std::smatch match;
            std::regex_search(header, match, start);
            absPos = std::stoi(match[1]); // read wig start position to absolute position
            
            if(argv[2] == NULL) {
                
                std::cout<<header<<"\n";
                
                uint8_t *values = new uint8_t[len*3];
                ifs.read(reinterpret_cast<char *>(&values), sizeof(uint8_t)*3*len);
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
                    
            } else if (argv[2] == std::string("--expand")){
                
                std::regex chrom("chrom\\=([^ ]*) ");
                std::regex_search(header, match, chrom);
                header = match[1];
                std::vector<uint8_t> kmerCov(k-1,0);
                std::vector<uint8_t> edgeCovFw(k-1,0);
                std::vector<uint8_t> edgeCovBw(k-1,0);
                
                for (uint64_t i = 0; i < len; ++i) { // loop through all position for this record
                    ifs.read(reinterpret_cast<char *>(&values), sizeof(uint8_t)*3);
                    
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
