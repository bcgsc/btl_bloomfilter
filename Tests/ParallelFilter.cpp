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
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace opt {
    unsigned kmerLen = 32;
    unsigned ibits = 8;
    unsigned nhash = 5;
}

using namespace std;

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
    for (size_t i = 0; i < seq.size() - opt::kmerLen + 1; i++) {
        string kmer = seq.substr(i,opt::kmerLen);
        myFilter.insert(kmer.c_str());
    }
        
    //string kmer = seq.substr(0,opt::kmerLen);
    //uint64_t fhVal, rhVal;
    //myFilter.insert(kmer.c_str(), fhVal, rhVal);
    //for (size_t i = 1; i < seq.size() - opt::kmerLen + 1; i++) {
    //    myFilter.insert(fhVal, rhVal, seq[i-1], seq[i+opt::kmerLen-1]);
    //}
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

void queryBf(const BloomFilter &myFilter, const char* faqFile) {
    size_t fHit=0;
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
            if(myFilter.contains(line.c_str())) {
                #pragma omp atomic
                ++fHit;
            }
            
                //__sync_add_and_fetch(&fHit, 1);
        }
    }
    uFile.close();
    cerr << "false hits = " << fHit << " " << setprecision(4) << fixed << (double)fHit/100000000.00 << "\n";
    //return fHit;
}

int main(int argc, const char* argv[]) {
    /*BloomFilter myFilter(40857600000, 2, 30);
    string mystr="AGAGACGTGCATCGGGTCATCAACCAATAT";
    myFilter.insert(mystr.c_str());
    if (myFilter.search(mystr.c_str()))
        cerr << mystr << " is in Bloom filter.\n";
    else
        cerr << mystr << " is not in Bloom filter.\n";
    return 0;*/
	if (argc <2) cerr << "error!\n";
    double sTime = omp_get_wtime();
    BloomFilter myFilter(800000000, opt::nhash, opt::kmerLen);
    loadBf(myFilter, argv[1]);
    cerr << "|popBF|=" << myFilter.getPop() << " ";
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    sTime = omp_get_wtime();
    queryBf(myFilter, argv[2]);
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    
    
    //myFilter.store("filter3.bf");
    
    /*BloomFilter filter2(40857600000, 5, 30, "filter1.bf");
    cerr << "|popBF|=" << filter2.getPop() << " ";
    cerr << setprecision(4) << fixed << omp_get_wtime() - sTime << "\n";
    
    filter2.store("filter2.bf");*/
    
    return 0;
}
