// New changes

#include <string>
#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <stdint.h>
#include "BloomFilterTest.hpp"
#include "BloomFilterInfo.hpp"
#include "seqgen.hpp"
#include <getopt.h>
//#include "city.h"
//#include "xxhash.h"

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
"Report bugs to hmohamadi@bcgsc.ca.\n";

namespace opt {
    unsigned threads;
    unsigned kmerLen;
    unsigned ibits = 8;
    unsigned nhash;
    size_t nquery;
    size_t squery;
    size_t ngene;
    size_t sgene;
    unsigned method;
    unsigned maxitr;
}

using namespace std;

static const char shortopts[] = "k:b:h:j:q:l:t:g:m:";

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
    { "max",	required_argument, NULL, 'm' },
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
    
    string kmer = seq.substr(0,opt::kmerLen);
    uint64_t fhVal, rhVal;
    myFilter.insert(kmer.c_str(), fhVal, rhVal);
    for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
        myFilter.insert(fhVal, rhVal, seq[i-1], seq[i+opt::kmerLen-1]);
    }
}

void loadSeqm(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        myFilter.insertMur(kmer.c_str());
    }
}

void loadSeqc(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        myFilter.insertCit(kmer.c_str());
    }
}

void loadSeqx(BloomFilter & myFilter, const string & seq) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        myFilter.insertXxh(kmer.c_str());
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
        if(good) {
            if(opt::method==0)
                loadSeqc(myFilter, line);
            else if(opt::method==1)
                loadSeqm(myFilter, line);
            else if(opt::method==2)
                loadSeqr(myFilter, line);
            else if(opt::method==3)
                loadSeqx(myFilter, line);
        }
    }
    uFile.close();
}

void querySeqr(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
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

void querySeqm(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(myFilter.containsMur(kmer.c_str())) {
            #pragma omp atomic
            ++fHit;
        }
    }
}

void querySeqc(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(myFilter.containsCit(kmer.c_str())) {
#pragma omp atomic
            ++fHit;
        }
    }
}

void querySeqx(BloomFilter & myFilter, const string & seq, size_t & fHit) {
    if (seq.size() < opt::kmerLen) return;
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i, opt::kmerLen);
        getCanon(kmer);
        if(myFilter.containsXxh(kmer.c_str())) {
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
            if(opt::method==0)
                querySeqc(myFilter, line, fHit);
            else if(opt::method==1)
                querySeqm(myFilter, line, fHit);
            else if(opt::method==2)
                querySeqr(myFilter, line, fHit);
            else if(opt::method==3)
                querySeqx(myFilter, line, fHit);
            
            #pragma omp atomic
            totKmer+=opt::squery-opt::kmerLen+1;
        }
    }
    uFile.close();
    cerr << "tkmer=" << totKmer << " ";
    cerr << "fhits=" << fHit << " %" << setprecision(4) << fixed << (double)fHit/(double)totKmer << " ";
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
            case 'm':
                arg >> opt::maxitr; break;
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
    
    
#ifdef _OPENMP
    omp_set_num_threads(opt::threads);
#endif

    //std::cerr<<"kmerl="<<opt::kmerLen<<"\n";
    //std::cerr<<"nhash="<<opt::nhash<<"\n";
    std::cerr<<"thread="<<opt::threads<<"\n";
    std::cerr<<"bit/i="<<opt::ibits<<"\n";

    std::cerr<<"nquery="<<opt::nquery<<"\n";
    std::cerr<<"squery="<<opt::squery<<"\n";
    std::cerr<<"ngene="<<opt::ngene<<"\n";
    std::cerr<<"sgene="<<opt::sgene<<"\n\n";

    
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
    string itm[4]={"city","murmur","rolling","xxhash"};
    
    for(unsigned trial=1; trial <= opt::maxitr; trial++) {
    makeGene(opt::ngene, opt::sgene);
    makeRead(opt::nquery, opt::squery);
    
    for (unsigned k=32; k<=160; k+=32) {
        opt::kmerLen = k;
        std::cerr<<"kmerl="<<opt::kmerLen<<"\n";
        for (unsigned i=1; i<6; i++) {
            opt::nhash = i;
            std::cerr<<"nhash="<<opt::nhash<<"\n";
            
            for(opt::method=0; opt::method<4; opt::method++) {
                std::cerr<<"method="<<itm[opt::method]<<"\n";
                double sTime = omp_get_wtime();
                BloomFilter myFilter(opt::ibits*opt::ngene*opt::sgene , opt::nhash, opt::kmerLen);
                loadBf(myFilter, geneName);
                cerr << "|popBF|=" << myFilter.getPop() << " ";
                cerr << "load_time=" <<setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
                
                sTime = omp_get_wtime();
                queryBf(myFilter, readName);
                cerr << "query_time=" << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
            }
            
            
            cerr << "nominal_fdr%=" << BloomFilterInfo::calcApproxFPR(opt::ibits*opt::ngene*opt::sgene, opt::ngene*opt::sgene,opt::nhash ) << "\n\n";
        }
    }
    }
    
    //myFilter.store("filter3.bf");
    
    /*BloomFilter filter2(40857600000, 5, 30, "filter1.bf");
    cerr << "|popBF|=" << filter2.getPop() << " ";
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    filter2.store("filter2.bf");*/
    
    return 0;
}
