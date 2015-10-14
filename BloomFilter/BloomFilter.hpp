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
#include "city.h"
#include <math.h>
#include <fstream>
#include <iostream>
#include <sys/stat.h>
#include <cstring>
#include <cassert>
#include <cstdlib>
#include <stdio.h>
#include <cstring>

using namespace std;

static const uint8_t bitsPerChar = 0x08;
static const unsigned char bitMask[0x08] = { 0x01, 0x02, 0x04, 0x08, 0x10, 0x20,
		0x40, 0x80 };

//TODO Work out better way to deal with kmerSize since conversion to kmerSize in bytes is needed

/*
 * For precomputing hash values. kmerSize is the number of bytes of the original string used.
 */
static inline vector<size_t> multiHash(const unsigned char* kmer, size_t num,
		unsigned kmerSize) {
	vector < size_t > tempHashValues(num);
	size_t kmerSizeInBytes = (kmerSize + 3) / 4;

	for (size_t i = 0; i < num; ++i) {
		tempHashValues[i] = CityHash64WithSeed(
				reinterpret_cast<const char*>(kmer), kmerSizeInBytes, i);
	}
	return tempHashValues;
}

class BloomFilter {
public:
	/* De novo filter constructor.
	 *
	 * preconditions:
	 * filterSize must be a multiple of 64
	 * kmerSize refers to the number of bases the kmer has
	 * k-mers supplied to this object should be binary (2 bits per base)
	 */
	BloomFilter(size_t filterSize, unsigned hashNum,
			unsigned kmerSize) :
			m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_kmerSizeInBytes(
					(kmerSize + 4 - 1) / 4) {
		initSize(m_size);
		memset(m_filter, 0, m_sizeInBytes);
	}

	/*
	 * Loads the filter (file is a .bf file) from path specified
	 */
	BloomFilter(size_t filterSize, unsigned hashNum,
			unsigned kmerSize, string const &filterFilePath) :
			m_size(filterSize), m_hashNum(hashNum), m_kmerSize(kmerSize), m_kmerSizeInBytes(
					(kmerSize + 4 - 1) / 4) {
		initSize(m_size);

		FILE *file = fopen(filterFilePath.c_str(), "rb");
		if (file == NULL) {
			cerr << "file \"" << filterFilePath << "\" could not be read."
					<< endl;
			exit(1);
		}

		long int lCurPos = ftell(file);
		fseek(file, 0, 2);
		size_t fileSize = ftell(file);
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

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	void insert(vector<size_t> const &precomputed) {

		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = precomputed.at(i) % m_size;
			__sync_or_and_fetch(&m_filter[normalizedValue / bitsPerChar],
					bitMask[normalizedValue % bitsPerChar]);
			//		m_filter[normalizedValue / bitsPerChar] |= bitMask[normalizedValue
			//				% bitsPerChar];
		}
	}

	void insert(const unsigned char* kmer) {
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = CityHash64WithSeed(
					reinterpret_cast<const char*>(kmer), m_kmerSizeInBytes, i)
					% m_size;
			__sync_or_and_fetch(&m_filter[normalizedValue / bitsPerChar],
					bitMask[normalizedValue % bitsPerChar]);
			//		m_filter[normalizedValue / bitsPerChar] |= bitMask[normalizedValue
			//				% bitsPerChar];
		}
	}

	/*
	 * Accepts a list of precomputed hash values. Faster than rehashing each time.
	 */
	bool contains(vector<size_t> const &values) const {
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = values.at(i) % m_size;
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
	bool contains(const unsigned char* kmer) const {
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t normalizedValue = CityHash64WithSeed(
					reinterpret_cast<const char*>(kmer), m_kmerSizeInBytes, i)
					% m_size;
			unsigned char bit = bitMask[normalizedValue % bitsPerChar];
			if ((m_filter[normalizedValue / bitsPerChar] & bit) != bit) {
				return false;
			}
		}
		return true;
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
		//write out each block
		myFile.write(reinterpret_cast<char*>(m_filter), m_sizeInBytes);

		myFile.close();
		assert(myFile);
	}

	unsigned getHashNum() const {
		return m_hashNum;
	}

	unsigned getKmerSize() const {
		return m_kmerSize;
	}

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

	uint8_t* m_filter;
	size_t m_size;
	size_t m_sizeInBytes;
	unsigned m_hashNum;
	unsigned m_kmerSize;
	unsigned m_kmerSizeInBytes;
};

#endif /* BLOOMFILTER_H_ */
