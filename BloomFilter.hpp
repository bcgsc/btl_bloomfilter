/*
 *
 * BloomFilter.hpp
 *
 *  Created on: Aug 10, 2012
 *      Author: cjustin
 */

#ifndef BLOOMFILTER_H_
#define BLOOMFILTER_H_
#include <string>
#include <vector>
#include <stdint.h>
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <cstring>
#include "rolling.h"

using namespace std;

struct FileHeader {
    char magic[8];
    uint32_t hlen;
    uint64_t size;
    uint32_t nhash;
    uint32_t kmer;
    double dFPR;
    double aFPR;
    double rFPR;
    uint64_t nEntry;
    uint64_t tEntry;
};


static const uint8_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
    0x40, 0x80 };

inline unsigned popCnt(unsigned char x) {
    return ((0x876543210 >>
             (((0x4332322132212110 >> ((x & 0xF) << 2)) & 0xF) << 2)) >>
            ((0x4332322132212110 >> (((x & 0xF0) >> 2)) & 0xF) << 2))
    & 0xf;
}

class BloomFilter {
public:

    /*
     * Default constructor.
    */
    BloomFilter() : m_filter(0), m_size(0), m_sizeInBytes(0), m_hashNum(0),
         m_kmerSize(0), m_dFPR(0), m_aFPR(0), m_rFPR(0), m_nEntry(0), m_tEntry(0) {}

    /* De novo filter constructor.
     *
     * preconditions:
     * filterSize must be a multiple of 64
     * kmerSize refers to the number of bases the kmer has
     */
    BloomFilter(size_t filterSize, unsigned hashNum, unsigned kmerSize) :
			m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_dFPR(
					0), m_aFPR(0), m_rFPR(0), m_nEntry(0), m_tEntry(0) {
        initSize(m_size);
        memset(m_filter, 0, m_sizeInBytes);
    }

    BloomFilter(const string &filterFilePath) {
        FILE *file = fopen(filterFilePath.c_str(), "rb");
	    if (file == NULL) {
		    cerr << "file \"" << filterFilePath << "\" could not be read."
					<< endl;
			exit(1);
		}

        loadHeader(file);

        long int lCurPos = ftell(file);
        fseek(file, 0, 2);
        size_t fileSize = ftell(file) - sizeof(struct FileHeader);
        fseek(file, lCurPos, 0);
        if (fileSize != m_sizeInBytes) {
            cerr << "Error: " << filterFilePath
            << " does not match size given by its information file. Size: "
            << fileSize << " vs " << m_sizeInBytes << " bytes." << endl;
            exit(1);
        }

        size_t countRead = fread(m_filter, fileSize, 1, file);
        if (countRead != 1 && fclose(file) != 0) {
            cerr << "file \"" << filterFilePath << "\" could not be read."
            << endl;
            exit(1);
        }
    }

    void loadHeader(FILE *file) {

        FileHeader header;
        fread(&header, sizeof(struct FileHeader), 1, file);
        char magic[9];
        strncpy(magic, header.magic, 8);
        magic[8] = '\0';

//        cerr << "Loading header... magic: " <<
//            magic << " hlen: " <<
//            header.hlen << " size: " <<
//            header.size << " nhash: " <<
//            header.nhash << " kmer: " <<
//            header.kmer << " dFPR: " <<
//            header.dFPR << " aFPR: " <<
//            header.aFPR << " rFPR: " <<
//            header.rFPR << " nEntry: " <<
//            header.nEntry << " tEntry: " <<
//            header.tEntry << endl;

        m_size = header.size;
        initSize(m_size);
        m_hashNum = header.nhash;
        m_kmerSize = header.kmer;
    }

    /*
     * Accepts a list of precomputed hash values. Faster than rehashing each time.
     */
    void insert(vector<size_t> const &precomputed) {

        //iterates through hashed values adding it to the filter
        for (size_t i = 0; i < m_hashNum; ++i) {
            size_t normalizedValue = precomputed.at(i) % m_size;
            __sync_or_and_fetch(&m_filter[normalizedValue / bitsPerChar],
                                bitMask[normalizedValue % bitsPerChar]);
        }
    }

    void insert(const char* kmer) {
        uint64_t hVal = getChval(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
            __sync_or_and_fetch(&m_filter[normalizedValue / bitsPerChar],
                                bitMask[normalizedValue % bitsPerChar]);
        }
    }

    /*
     * Accepts a list of precomputed hash values. Faster than rehashing each time.
     */
    bool contains(vector<size_t> const &precomputed) const {
        for (size_t i = 0; i < m_hashNum; ++i) {
            size_t normalizedValue = precomputed.at(i) % m_size;
            unsigned char bit = bitMask[normalizedValue % bitsPerChar];
            if ((m_filter[normalizedValue / bitsPerChar] & bit) != bit) {
                return false;
            }
        }
        return true;
    }

    /*
     * Single pass filtering, computes hash values on the fly
     */
    bool contains(const char* kmer) const {
        uint64_t hVal = getChval(kmer, m_kmerSize);
        for (unsigned i = 0; i < m_hashNum; i++) {
            size_t normalizedValue = (rol(varSeed, i) ^ hVal) % m_size;
            unsigned char bit = bitMask[normalizedValue % bitsPerChar];
            if ((m_filter[normalizedValue / bitsPerChar] & bit) == 0)
                return false;
        }
        return true;
    }

    void writeHeader(ofstream &out) const {
        FileHeader header;
        strncpy(header.magic, "BlOOMFXX", 8);
        char magic[9];
        strncpy(magic, header.magic, 8);
        magic[8] = '\0';

        header.hlen = sizeof(struct FileHeader);
        header.size = m_size;
        header.nhash = m_hashNum;
        header.kmer = m_kmerSize;
        header.dFPR = m_dFPR;
        header.aFPR = m_aFPR;
        header.rFPR = m_rFPR;
        header.nEntry = m_nEntry;
        header.tEntry = m_tEntry;

//        cerr << "Writing header... magic: "
//            << magic << " hlen: "
//            << header.hlen << " size: "
//            << header.size << " nhash: "
//            << header.nhash << " kmer: "
//            << header.kmer << " dFPR: "
//            << header.dFPR << " aFPR: "
//            << header.aFPR << " rFPR: "
//            << header.rFPR << " nEntry: "
//            << header.nEntry << " tEntry: "
//            << header.tEntry << endl;

        out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));
    }


    /*
     * Stores the filter as a binary file to the path specified
     * Stores uncompressed because the random data tends to
     * compress poorly anyway
     */
    void storeFilter(string const &filterFilePath) const {
        ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

        cerr << "Storing filter. Filter is " << m_sizeInBytes << "bytes."
        << endl;

        assert(myFile);
        writeHeader(myFile);

        //write out each block
        myFile.write(reinterpret_cast<char*>(m_filter), m_sizeInBytes);

        myFile.close();
        assert(myFile);
    }

    size_t getPop() const {
        size_t i, popBF=0;
#pragma omp parallel for reduction(+:popBF)
        for(i=0; i<(m_size + 7)/8; i++)
            popBF = popBF + popCnt(m_filter[i]);
        return popBF;
    }

    unsigned getHashNum() const {
        return m_hashNum;
    }

    unsigned getKmerSize() const {
        return m_kmerSize;
    }

    void setdFPR(double value) {
        m_dFPR = value;
    }

    void setaFPR(double value) {
        m_aFPR = value;
    }

    void setrFPR(double value) {
        m_rFPR = value;
    }

    void setnEntry(uint64_t value) {
        m_nEntry = value;
    }
    
    void settEntry(uint64_t value) {
        m_tEntry = value;
    }

    size_t getFilterSize() const { return m_size; }

    ~BloomFilter() {
        delete[] m_filter;
    }
private:
    BloomFilter(const BloomFilter& that); //to prevent copy construction

    /*
     * Checks filter size and initializes filter
     */
    void initSize(size_t size) {
        if (size % 8 != 0) {
            cerr << "ERROR: Filter Size \"" << size
            << "\" is not a multiple of 8." << endl;
            exit(1);
        }
        m_sizeInBytes = size / bitsPerChar;
        m_filter = new unsigned char[m_sizeInBytes];
    }

    /*
	 * Only returns multiples of 64 for filter building purposes
	 * Is an estimated size using approximations of FPR formula
	 * given the number of hash functions
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	size_t calcOptimalSize(size_t entries, float fpr, unsigned hashNum) const {
		size_t non64ApproxVal = size_t(
				-double(entries) * double(hashNum)
						/ log(1.0 - pow(fpr, float(1 / (float(hashNum))))));

		return non64ApproxVal + (64 - non64ApproxVal % 64);
	}

    uint8_t* m_filter;
    size_t m_size;
    size_t m_sizeInBytes;
    unsigned m_hashNum;
    unsigned m_kmerSize;
    double m_dFPR;
    double m_aFPR;
    double m_rFPR;
    uint64_t m_nEntry;
    uint64_t m_tEntry;
};

#endif /* BLOOMFILTER_H_ */
