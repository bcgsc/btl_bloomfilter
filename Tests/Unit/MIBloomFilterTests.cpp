/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN
/* lightweight unit test framework */
#include "catch.hpp"
#include <sdsl/int_vector.hpp>
#include <stdio.h>

using namespace std;

#include "MIBloomFilter.hpp"
#include "ntHashIterator.hpp"

TEST_CASE("Test for bit vector acting as bloom filter", "[MIBloomFilter]")
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
    unsigned estEntries = seq1.length() + seq2.length() + seq3.length() - 3*k + 3;

    SECTION("1 Hash Function")
    {
        const unsigned numHashes = 1;
        size_t bvSize = MIBloomFilter<size_t>::calcOptimalSize(estEntries, numHashes, 0.5);    
        sdsl::bit_vector bv(bvSize);
        REQUIRE((bv.size() > 0));
        
        for (unsigned i = 0; i < 3; i++){ 
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                collisionCount += MIBloomFilter<size_t>::insert(bv,(uint64_t*) *itr, numHashes);
                bitCount++;
                itr.next();
            }
        }
        for(size_t i = 0; i < bvSize; i++)
            bitSum += bv[i];
        REQUIRE(collisionCount == numHashes * bitCount - bitSum);
    }

    SECTION("5 Hash Functions")
    {
        const unsigned numHashes = 5;
        size_t bvSize = MIBloomFilter<size_t>::calcOptimalSize(estEntries, numHashes, 0.5);    
        sdsl::bit_vector bv(bvSize);
        REQUIRE((bv.size() > 0));

        for (unsigned i = 0; i < 3; i++){ 
            ntHashIterator itr(seqArr[i], numHashes, k);
            while(itr != itr.end()){
                collisionCount += MIBloomFilter<size_t>::insert(bv,(uint64_t*) *itr, numHashes);
                bitCount++;
                itr.next();
            }
        }
        for(size_t i = 0; i < bvSize; i++)
            bitSum += bv[i];
        REQUIRE(collisionCount == numHashes * bitCount - bitSum);
    }
}


