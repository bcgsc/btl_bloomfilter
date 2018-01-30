/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN
/* lightweight unit test framework */
#include "catch.hpp"
#include <sdsl/int_vector.hpp>
#include <stdio.h>

using namespace std;

#include "MIBloomFilter.hpp"
#include "ntHashIterator.hpp"
#include "MIBloomFilterUtil.hpp"

TEST_CASE("Test loading", "[MIBloomFilter]")
{
    /* Setup */
    const string seq1= "TAACGGGCGATTCTATAAGATTGCACATTGCGTCTACTTATAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTGCCTACTATCCTTAAACGCATATCT";
    const string seq2= "CTCCGGCAAGCAATTATGAACAACGCAAGGATCGGCGATATAAACAGAGAAACGGCTGATTACACTTGTTCGTGTGGTATCGCTAAATAGCCTCGCGGAG";
    const string seq3= "GTCGGCCCCATCAGTAGCCCGAATATGTCGCTTTACGGGTCCTGGGCCGGGGTGCGATACCTTGCAGAAATCGAGGCCGTTCGTTAATTCCTGTTGCATT";
    const string seqArr[3] = {seq1, seq2, seq3};

    unsigned bitCount = 0;
    unsigned collisionCount = 0;
    unsigned bitSum = 0;
    const unsigned k = 15;
    size_t estEntries = seq1.length() + seq2.length() + seq3.length() - 3*k + 3;

    SECTION("Bit vector for 1 hash")
    {   //construct a bit vector
        const unsigned numHashes = 1;
        size_t bvSize = MIBloomFilter<size_t>::calcOptimalSize(estEntries, numHashes, 0.5);    
        sdsl::bit_vector bv(bvSize);
        REQUIRE((bv.size() > 0));

        //Using the static function insert the values into the bit vector
        for (unsigned i = 0; i < 3; ++i){ 
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                collisionCount += MIBloomFilter<size_t>::insert(bv,(uint64_t*) *itr, numHashes);
                ++bitCount;
                ++itr;
            }
        }

        //Calulate the total number of i in the bit vector
        for(size_t i = 0; i < bvSize; ++i)
            bitSum += bv[i];
        //Determines if the non collided inserts are in the vector
        REQUIRE(collisionCount == numHashes * bitCount - bitSum);
    }

    //Virtually identical to the above section, but uses 5 hashes this time
    SECTION("Bit vector for 5 hashes")
    {
        const unsigned numHashes = 5;
        size_t bvSize = MIBloomFilter<size_t>::calcOptimalSize(estEntries, numHashes, 0.5);    
        sdsl::bit_vector bv(bvSize);
        REQUIRE((bv.size() > 0));

        for (unsigned i = 0; i < 3; ++i){ 
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                collisionCount += MIBloomFilter<size_t>::insert(bv,(uint64_t*) *itr, numHashes);
                ++bitCount;
                ++itr;
            }
        }
        for(size_t i = 0; i < bvSize; ++i)
            bitSum += bv[i];
        REQUIRE(collisionCount == numHashes * bitCount - bitSum);
    }

    SECTION("Bloom Filter for 1 Hash")
    {   //Construct bit vector for filter
        const unsigned numHashes = 1;
        bool saturated = true;
        bool contains = false;
        size_t bvSize = MIBloomFilter<size_t>::calcOptimalSize(estEntries, numHashes, 0.5);
        sdsl::bit_vector bv(bvSize);
        REQUIRE(bv.size() > 0);
        for (unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                MIBloomFilter<size_t>::insert(bv,(uint64_t*) *itr, numHashes);
                ++itr;
            }
        }
        //Constructs a filter and loads all values into filter
        MIBloomFilter<uint8_t> filter(numHashes, k, bv);
        unsigned max = 1;
        for (unsigned i = 0; i < 3; ++i){
            unsigned ID = i + 1;
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                saturated = true;
                filter.insert(*itr, ID, max, saturated);
                ++itr;
            }
            vector<uint8_t> results = filter.at(*itr, saturated);
        }

        //Determines if ID is present in the filter if filter is not saturated
        for (unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                saturated = true;
                contains = false;
                vector<uint8_t> results = filter.at(*itr, saturated);
                if(!saturated){
                    for(int ID = 1; ID <= 3; ++ID){
                        if(results[0] == ID){
                            contains = true;
                            break;
                        }
                    }
                    REQUIRE(contains);
                }
                ++itr;
            }
        }
        REQUIRE(filter.getPopNonZero()!=0);
    }
    SECTION("Bloom Filter for 5 Hashes")
    {
        const unsigned numHashes = 5;
        bool saturated = true;
        size_t bvSize = MIBloomFilter<size_t>::calcOptimalSize(estEntries, numHashes, 0.5);
        sdsl::bit_vector bv(bvSize);
        REQUIRE(bv.size() > 0);
        for (unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                MIBloomFilter<size_t>::insert(bv,(uint64_t*) *itr, numHashes);
                ++itr;
            }
        }
        //Constructs a filter and loads all values into filter
        MIBloomFilter<uint8_t> filter(numHashes, k, bv);
        for(unsigned max = 1; max <= numHashes; ++max){
            for (unsigned i = 0; i < 3; ++i){
                unsigned ID = i + 1;
                ntHashIterator itr(seqArr[i], numHashes, k);
                while(itr != itr.end()){
                    saturated = true;
                    filter.insert(*itr, ID, max, saturated);
                    ++itr;
                }
                vector<uint8_t> results = filter.at(*itr, saturated);
            }
        }
        //Determines if ID is present in the filter if filter is not saturated
        for (unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                saturated = true;
                vector<uint8_t> results = filter.at(*itr, saturated);
                if(!saturated){
                    int presenceCount = 0; //Number of times the ID appears in the vector
                    for(unsigned resI = 0; resI < results.size(); ++resI)
                        for(int ID = 1; ID <= 3; ++ID)
                            if(results[resI] == ID)
                                ++presenceCount;
                    REQUIRE(presenceCount == results.size());
                }
                ++itr;
            }
        }
        REQUIRE(filter.getPopNonZero()!=0);
    }

}

TEST_CASE("Testing Query Function", "[MIBloomFilter]")
{
    /* Setup */
    const string seq1= "TAACGGGCGATTCTATAAGATTGCACATTGCGTCTACTTATAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTGCCTACTATCCTTAAACGCATATCT";
    const string seq2= "CTCCGGCAAGCAATTATGAACAACGCAAGGATCGGCGATATAAACAGAGAAACGGCTGATTACACTTGTTCGTGTGGTATCGCTAAATAGCCTCGCGGAG";
    const string seq3= "GTCGGCCCCATCAGTAGCCCGAATATGTCGCTTTACGGGTCCTGGGCCGGGGTGCGATACCTTGCAGAAATCGAGGCCGTTCGTTAATTCCTGTTGCATT";
    const string seqArr[3] = {seq1, seq2, seq3};

    unsigned bitCount = 0;
    unsigned collisionCount = 0;
    unsigned bitSum = 0;
    const unsigned k = 25;
    size_t estEntries = seq1.length() + seq2.length() + seq3.length() - 3*k + 3;

    const unsigned numHashes = 5;
    bool saturated = true;
    bool contains = false;
    size_t bvSize = MIBloomFilter<size_t>::calcOptimalSize(estEntries, numHashes, 0.5);
    sdsl::bit_vector bv(bvSize);
    REQUIRE(bv.size() > 0);
    for (unsigned i = 0; i < 3; ++i){
        ntHashIterator itr(seqArr[i], numHashes, k);
        while(itr != itr.end()){
            MIBloomFilter<size_t>::insert(bv,(uint64_t*) *itr, numHashes);
            ++itr;
        }
    }
    //Constructs a filter and loads all values into filter
	MIBloomFilter<uint8_t> filter(numHashes, k, bv);
	for (unsigned max = 1; max <= numHashes; ++max) {
		for (unsigned i = 0; i < 3; ++i) {
			unsigned ID = i + 1;
			ntHashIterator itr(seqArr[i], numHashes, k);
			while (itr != itr.end()) {
				saturated = true;
				filter.insert(*itr, ID, max, saturated);
				++itr;
			}
		}
	}
    
    const vector<double> frameProb = MIBloomFilterUtil::calcPerFrameProb(filter, (unsigned char) 3);
    
    SECTION("Testing Identical Sequences")
    {           
        for(unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(seqArr[i], numHashes, k);
            unsigned ID = i + 1;
            vector<uint8_t> results = MIBloomFilterUtil::query(filter, itr, frameProb);
            REQUIRE(results.size() == 1);
            REQUIRE(results[0] == ID);
        }
    }

    SECTION("Testing Subsequences")
    {   //Subsequences of seq1
        const string sub1no1 = "TAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTGCCTACTATCCTTAAACGCATATCT";
        const string sub1no2 = "TAACGGGCGATTCTATAAGATTGCACATTGCGTCTACTTATAAGATGTC";
        const string sub1no3 = "TAAGATGTCTCAACGGCATGCGCAACTTGTGAAGTGCCTACTATCC";
        //Subsequences of seq2
        const string sub2no1= "TCCGGCAAGCAATTATGAACAACGCAAGGATCGGCGATATAAACAGAGAAACGGCTGATTACACTTGTTCGTGTGGTATCGCTAAATAGCCTCGCGGA";
        const string sub2no2= "ACAACGCAAGGATCGGCGATATAAACAGAGAAACGGCTGATTACACTTGTTCGTGTGGTATCGCTAAATAGCCTCGCGGAG";
        const string sub2no3= "CTCCGGCAAGCAATTATGAACAACGCAAGGATCGGCGATATAAACAGAGAAACGGCTGAT";
        //Subsequences of seq3
        const string sub3no1 = "TTACGGGTCCTGGGCCGGGGTGCGATACCTTGCAGAAATCGAGGCCGTTCGTTAATTCCTGTTGCATT";
        const string sub3no2 = "GTCGGCCCCATCAGTAGCCCGAATATGTCGCTTTACGGGTCCTGGGCCGGGGTGCGATACCTTGCAGAAATCGAGGC";
        const string sub3no3 = "CGGGTCCTGGGCCGGGGTGCGATACCTTGCAGAAATCGAGGCCGTTCGTT";
        const string subArr[9] = {sub1no1, sub1no2, sub1no3, sub2no1, sub2no2, sub2no3, sub3no1, sub3no2, sub3no3};
        //Checks if sub1noX can be found in the filter
        for(unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(subArr[i], numHashes, k);
            vector<uint8_t> results = MIBloomFilterUtil::query(filter, itr, frameProb);
            REQUIRE(results.size() == 1);
            REQUIRE(results[0] == 1);
        }

        //Checks if sub2noX can be found in the filter
        for(unsigned i = 3; i < 6; ++i){
            ntHashIterator itr(subArr[i], numHashes, k);
            vector<uint8_t> results = MIBloomFilterUtil::query(filter, itr, frameProb);
            REQUIRE(results.size() == 1);
            REQUIRE(results[0] == 2);
        }

        //Checks if sub3noX can be found in the filter
        for(unsigned i = 6; i < 9; ++i){
            ntHashIterator itr(subArr[i], numHashes, k);
            vector<uint8_t> results = MIBloomFilterUtil::query(filter, itr, frameProb);
            REQUIRE(results.size() == 1);
            REQUIRE(results[0] == 3);
        }
    }

    SECTION("Testing Sequences not in Filter")
    {
        const string badSeq1 = "ATGTTCACCTATCTACTACCCATCCCCGGAGGTTAAGTAGGTTGTGAGATGCGGGAGAGGTTCTCGATCTTCCCG";
        const string badSeq2 = "TAGAGCGGGGCTGTTGACGTTTGGAGTTGAAAAAATCTAATATTCCAATC";
        const string badSeq3 = "GGCTTCAACGTGCACCACCGCAGGCGGCTGACGAGGGGCTCACACCGAGAAAGTAGACTGTTGCGCGTTGGGGGT";
        const string badArr[3] = {badSeq1, badSeq2, badSeq3};
        for(unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(badArr[i], numHashes, k);
            vector<uint8_t> results = MIBloomFilterUtil::query(filter, itr, frameProb);
            REQUIRE(results.size() == 0);
        }
    }

}
