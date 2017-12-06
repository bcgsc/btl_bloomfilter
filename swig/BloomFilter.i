%module BloomFilter
%include "std_string.i"
%include "stdint.i"
%include "std_vector.i"
/*%apply const uint64_t& { uint64_t & };*/
/*%apply uint64_t& INOUT {uint64_t& fhVal, uint64_t& rhVal};*/
namespace std {
   %template(SizetVector) vector<size_t>;
}

%{
#include "../KmerBloomFilter.hpp"
#include "../ntHashIterator.hpp"
#include "../BloomFilterUtil.h"
%}

%rename(BloomFilter) KmerBloomFilter;

using namespace std;

class KmerBloomFilter {
public:
        KmerBloomFilter();
        ~KmerBloomFilter();
        KmerBloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize);
        KmerBloomFilter(const string &filterFilePath);

        void insert(vector<size_t> const &precomputed);
        void insert(const char* kmer);

        bool contains(vector<size_t> const &values);
        bool contains(const char* kmer);

        void storeFilter(string const &filterFilePath);
        size_t getPop();
        unsigned getHashNum();
        unsigned getKmerSize();
        size_t getFilterSize();
};

/*
class ntHashIterator {
public:
    ntHashIterator();
    ~ntHashIterator();
    ntHashIterator(const string& seq, unsigned numHashes, unsigned k);
 
    const uint64_t* operator*();
    const uint64_t* operator->();

    bool operator==(const ntHashIterator& it);
    bool operator!=(const ntHashIterator& it);
    
    ntHashIterator& operator++();
    static const ntHashIterator end();
};
*/

void insertSeq(KmerBloomFilter &bloom, const string& seq, unsigned numHashes, unsigned k);
