/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN
/* lightweight unit test framework */
#include "catch.hpp"
#include <sdsl/int_vector.hpp>
#include <stdio.h>

using namespace std;

#include "MIBloomFilter.hpp"
#include "ntHashIterator.hpp"

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
    const unsigned k = 4;
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
                if(!filter.insert(*itr, ID, max, saturated))
                    ++collisionCount;
                else
                    ++bitCount;
                ++itr;
            }
            vector<uint8_t> vec = filter.at(*itr, saturated);
        }

        //Determines if ID is present in the filter if filter is not saturated
        for (unsigned i = 0; i < 3; ++i){
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                saturated = true;
                contains = false;
                vector<uint8_t> vec = filter.at(*itr, saturated);
                if(!saturated){
                    for(int ID = 1; ID <= 3; ++ID){
                        if(vec[0] == ID){
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
}
