// New changes

#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdint.h>
#include "BloomFilter.hpp"
#include "seqgen.hpp"
#include <getopt.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM "ParallelFilter"

static const char VERSION_MESSAGE[] =
PROGRAM " Version 1.0.0 \n"
"Written by Hamid Mohamadi.\n"
"Copyright 2015 Canada's Michael Smith Genome Science Centre\n";

static const char USAGE_MESSAGE[] =
"Usage: " PROGRAM " [OPTION]... QUERY\n"
"Dispatch the sequences of the files QUERY based on the Bloom filter of the file TARGET.\n"
"\n"
" Options:\n"
"\n"
"  -p, --partition=N       divide reference to N partitions\n"
"  -j, --threads=N         use N parallel threads [partitions]\n"
"  -l, --alen=N            the minimum alignment length [20]\n"
"  -b, --bmer=N            size of a bmer [3*alen/4]\n"
"  -s, --step=N            step size used when breaking a query sequence into bmers [bmer]\n"
"  -h, --hash=N            use N hash functions for Bloom filter [6]\n"
"  -i, --bit=N             use N bits for each item in Bloom filter [8]\n"
"      --se                single-end library\n"
"      --fq                dispatch reads in fastq format\n"
"      --help              display this help and exit\n"
"      --version           output version information and exit\n"
"\n"
"Report bugs to hmohamadi@bcgsc.ca.\n";

namespace opt {
    unsigned threads;
    unsigned kmerLen=64;
    unsigned ibits = 8;
    unsigned nhash = 5;
    size_t nquery=100000000;
    size_t squery=200;
    size_t ngene=100;
    size_t sgene=5000000;
}

using namespace std;

static const char shortopts[] = "k:b:h:j:q:l:t:g:";

enum { OPT_HELP = 1, OPT_VERSION };

static const struct option longopts[] = {
    { "threads",	required_argument, NULL, 'j' },
    { "kmer",	required_argument, NULL, 'k' },
    { "qnum",	required_argument, NULL, 'q' },
    { "qlen",	required_argument, NULL, 'l' },
    { "bit",	required_argument, NULL, 'b' },
    { "hash",	required_argument, NULL, 'h' },
    { "tnum",	required_argument, NULL, 't' },
    { "tlen",	required_argument, NULL, 'g' },
    { "help",	no_argument, NULL, OPT_HELP },
    { "version",	no_argument, NULL, OPT_VERSION },
    { NULL, 0, NULL, 0 }
};

static const unsigned char b2r[256] = {
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //0
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //1
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //2
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //3
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //4   'A' 'C' 'G'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //5   'T'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'T', 'N', 'G', 'N', 'N', 'N', 'C', //6   'a' 'c' 'g'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'A', 'N', 'N', 'N', //7   't'
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //8
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //9
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //10
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //11
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //12
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //13
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //14
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N',
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N', //15
    'N', 'N', 'N', 'N', 'N', 'N', 'N', 'N'
};

void getCanon(std::string &bMer) {
    int p=0, hLen=(opt::kmerLen-1)/2;
    while (bMer[p] == b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        ++p;
        if(p>=hLen) break;
    }
    if (bMer[p] > b2r[(unsigned char)bMer[opt::kmerLen-1-p]]) {
        for (int lIndex = p, rIndex = opt::kmerLen-1-p; lIndex<=rIndex; ++lIndex,--rIndex) {
            char tmp = b2r[(unsigned char)bMer[rIndex]];
            bMer[rIndex] = b2r[(unsigned char)bMer[lIndex]];
            bMer[lIndex] = tmp;
        }
    }
}

void loadSeq(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        myFilter.insert(kmer.c_str());
    }
}

void loadSeqr(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    
    /*for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i,opt::kmerLen);
        myFilter.insert(kmer.c_str());
    }*/
        
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal;
    myFilter.insert(kmer.c_str(), fhVal, rhVal);
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        myFilter.insert(fhVal, rhVal, seq[i-1], seq[i+opt::kmerLen-1]);
    }
}


void loadBf(BloomFilter &myFilter, const char* faqFile) {
    ifstream uFile(faqFile);
    bool good = true;
    #pragma omp parallel
    for(string line, hline; good;) {
        #pragma omp critical(uFile)
        {
            good = getline(uFile, hline);
            good = getline(uFile, line);
            //good = getline(uFile, hline);
            //good = getline(uFile, hline);
        }
        if(good) loadSeqr(myFilter, line);
    }
    uFile.close();
}

void querySeqr(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    
    /*for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
     string kmer = seq.substr(i,opt::kmerLen);
     myFilter.insert(kmer.c_str());
     }*/
    
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal;
    if(myFilter.contains(kmer.c_str(), fhVal, rhVal)) {
        #pragma omp atomic
        ++fHit;
    }
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        if(myFilter.contains(fhVal, rhVal, seq[i-1], seq[i+opt::kmerLen-1])) {
            #pragma omp atomic
            ++fHit;
        }
            
    }
}

void queryBf(BloomFilter &myFilter, const char* faqFile) {
    size_t fHit=0,totKmer=0;
    ifstream uFile(faqFile);
    bool good = true;
#pragma omp parallel
    for(string line, hline; good;) {
#pragma omp critical(uFile)
        {
            good = getline(uFile, hline);
            good = getline(uFile, line);
            //good = getline(uFile, hline);
            //good = getline(uFile, hline);
        }
        if(good) {
            //getCanon(line);
            /*if(myFilter.contains(line.c_str())) {
                #pragma omp atomic
                ++fHit;
            }*/
            querySeqr(myFilter, line, fHit);
            #pragma omp atomic
            totKmer+=opt::squery-opt::kmerLen+1;
            
                //__sync_add_and_fetch(&fHit, 1);
        }
    }
    uFile.close();
    cerr << "totKmer = " << totKmer << "\n";
    cerr << "false hits = " << fHit << " " << setprecision(4) << fixed << (double)fHit/(double)totKmer << "\n";
    //return fHit;
}


int main(int argc, char** argv) {

    bool die = false;
    for (int c; (c = getopt_long(argc, argv, shortopts, longopts, NULL)) != -1;) {
        std::istringstream arg(optarg != NULL ? optarg : "");
        switch (c) {
            case '?':
                die = true; break;
            case 'j':
                arg >> opt::threads; break;
            case 'b':
                arg >> opt::ibits; break;
            case 'q':
                arg >> opt::nquery; break;
            case 'l':
                arg >> opt::squery; break;
            case 't':
                arg >> opt::ngene; break;
            case 'g':
                arg >> opt::sgene; break;
            case 'h':
                arg >> opt::nhash; break;
            case 'k':
                arg >> opt::kmerLen; break;
            case OPT_HELP:
                std::cerr << USAGE_MESSAGE;
                exit(EXIT_SUCCESS);
            case OPT_VERSION:
                std::cerr << VERSION_MESSAGE;
                exit(EXIT_SUCCESS);
        }
        if (optarg != NULL && !arg.eof()) {
            std::cerr << PROGRAM ": invalid option: `-"
            << (char)c << optarg << "'\n";
            exit(EXIT_FAILURE);
        }
    }
    //cerr << "argc=" << argc << "\noptind=" << optind << "\n"; exit(0);
    if (argc - optind != 2) {
        std::cerr << PROGRAM ": missing arguments\n";
        die = true;
    }
    
    if (die) {
        std::cerr << "Try `" << PROGRAM
        << " --help' for more information.\n";
        exit(EXIT_FAILURE);
    }
    

    std::cerr<<"kmerl="<<opt::kmerLen<<"\n";
    std::cerr<<"nhash="<<opt::nhash<<"\n";
    std::cerr<<"bit/i="<<opt::ibits<<"\n";

    std::cerr<<"nquery="<<opt::nquery<<"\n";
    std::cerr<<"squery="<<opt::squery<<"\n";
    std::cerr<<"ngene="<<opt::ngene<<"\n";
    std::cerr<<"sgene="<<opt::sgene<<"\n";

    
    const char *geneName(argv[argc-2]);
    const char *readName(argv[argc-1]);

    
    /*BloomFilter myFilter(40857600000, 2, 30);
    string mystr="AGAGACGTGCATCGGGTCATCAACCAATAT";
    myFilter.insert(mystr.c_str());
    if (myFilter.search(mystr.c_str()))
        cerr << mystr << " is in Bloom filter.\n";
    else
        cerr << mystr << " is not in Bloom filter.\n";
    return 0;*/
    makeGene(opt::ngene, opt::sgene+opt::kmerLen-1);
    makeRead(opt::nquery, opt::squery);
    
    double sTime = omp_get_wtime();
    BloomFilter myFilter(opt::ibits*opt::ngene*opt::sgene, opt::nhash, opt::kmerLen);
    loadBf(myFilter, geneName);
    cerr << "|popBF|=" << myFilter.getPop() << " ";
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    sTime = omp_get_wtime();
    queryBf(myFilter, readName);
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    
    
    //myFilter.store("filter3.bf");
    
    /*BloomFilter filter2(40857600000, 5, 30, "filter1.bf");
    cerr << "|popBF|=" << filter2.getPop() << " ";
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    filter2.store("filter2.bf");*/
    
    return 0;
}
