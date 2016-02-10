#ifndef BLOOM_FILTER_UTIL
#define BLOOM_FILTER_UTIL 1

#include "RollingHash.h"
#include "BloomFilter.hpp"
#include "RollingHashIterator.h"

using namespace std;

void insertSeq(BloomFilter &bloom, const string &seq, unsigned hashNum, unsigned kmerSize) {
    RollingHashIterator itr(seq, hashNum, kmerSize); 
    while (itr != itr.end()) {
        bloom.insert(*itr);
        itr++;
    }
}

#endif
