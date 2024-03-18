#include <stdlib.h>
#include <map>
#include <cstdio>

#include "validate.h"

int main(int, char **argv) {
    std::cout << "WARNING: only run this program if gfastats is in a working state" << std::endl;
    std::cout << "WARNING: previous validate files will be deleted" << std::endl;
    std::cout << "continue? (Y/N) ";
    std::string input;
    std::cin >> input;
    if(input != "Y" && input != "y") {
        std::cout << "validate generation cancelled" << std::endl;
        std::exit(0);
    }
    std::cout << "deleting old validate files..." << std::endl;

    for(auto &file : list_dir("validateFiles")) {
        if(getFileExt(file) != "tst") continue; // dont delete README
        file = "validateFiles/"+file;
        if(remove(file.c_str()) != 0) {
            std::cerr << "error deleting <" << file << ">" << std::endl;
            return -1;
        }
    }

    std::cout << "generating new validate files..." << std::endl;

    std::string exePath = getExePath(argv[0]);

    const std::map<std::set<std::string>, std::vector<std::string>> ext_args = {
        {{"fasta", "fasta.gz", "fastq", "fastq.gz", "gfa"}, {"-r testFiles/random1.fastq", "-r testFiles/random2.fastq", "-r testFiles/random1.fastq.gz", "-r testFiles/random1.fastq testFiles/random2.fastq", "-r testFiles/random1.fastq.gz testFiles/random2.fastq.gz"}}
    //  {{set of test file extensions}, {list of command line args to run with}}
    };
    
    const std::set<std::string> excludeExt {};
    const std::set<std::string> excludeFile {"random4.fasta", "random4.fastq", "to_correct.fasta", "to_correct.fastq", "decompressor1.fasta"};

    std::map<std::set<std::string>, std::vector<std::string>> file_args = {
        {{"random1.fasta"}, {"-r testFiles/random3.N.fastq", "-d testFiles/test1.kreeq", "-d testFiles/test2.kreeq"}},
        {{"random4.fasta"}, {"-r testFiles/random4.fastq -k3"}},
        {{"to_correct.fasta"}, {"-r testFiles/to_correct.fastq", "-r testFiles/to_correct.fastq -o gfa", "-r testFiles/to_correct.fastq -o vcf"}}
    //  {{set of test file paths}, {list of command line args to run with}}
    };

    for(const std::string &file : list_dir("testFiles")) {
        std::string ext = getFileExt(file);
        if(excludeExt.count(ext)) continue;
        if(excludeFile.count(file)) continue;
        for(auto pair : ext_args) {
            if(!pair.first.count(ext)) continue;
            for(auto args : pair.second) {
                genTest(exePath, file, args, "validate");
            }
        }
    }

    std::fstream fstream;
    for(const auto &pair : file_args) {
        for(const std::string &file : pair.first) {
            fstream.open("testFiles/"+file);
            if(!fstream) continue;
            fstream.close();
            for(const std::string &args : pair.second) {
                genTest(exePath, file, args, "validate");
            }
        }
    }
    
    file_args = {
        {{"-d testFiles/test1.kreeq testFiles/test2.kreeq"}, {""}}
    //  {{set of test file paths}, {list of command line args to run with}}
    };
    
    for(const auto &pair : file_args) {
        for(const std::string &file : pair.first) {
            for(const std::string &args : pair.second) {
                genTestUnion(exePath, file, args, "union");
            }
        }
    }

    std::exit(EXIT_SUCCESS);
}
