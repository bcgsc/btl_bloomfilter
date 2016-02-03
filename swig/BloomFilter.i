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
#include "../BloomFilter.hpp"
#include "../RollingHashIterator.h"
%}


using namespace std;

class BloomFilter {
public:
        BloomFilter();
        ~BloomFilter();
        BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize);
        BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize, const string &filterFilePath);

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

class RollingHashIterator {
public:
    RollingHashIterator();
    ~RollingHashIterator();
    RollingHashIterator(const string& seq, unsigned numHashes, unsigned k);
 
    const std::vector<size_t>& operator*();
    const std::vector<size_t>* operator->();

    bool operator==(const RollingHashIterator& it);
    bool operator!=(const RollingHashIterator& it);
    vector<size_t> * getNext();
/*  RollingHashIterator& operator++();
    RollingHashIterator operator++(int);
    static const RollingHashIterator end();*/
    unsigned pos();
    string kmer();
};
