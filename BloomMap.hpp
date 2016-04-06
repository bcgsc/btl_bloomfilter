/*
 * BloomMap.hpp
 *
 *  Created on: Dec 17, 2015
 *      Author: gjahesh
 */

#ifndef BLOOMMAP_HPP_
#define BLOOMMAP_HPP_

#include <vector>
#include <omp.h>
#include <fstream>
#include <iostream>
#include <stdio.h>

static const uint8_t bitsPerChar = 0x08;
template<typename T>

class BloomMap {

public:

	struct FileHeader {
		char magic[8];
		uint32_t hlen;
		uint64_t size;
		uint32_t nhash;
		uint32_t kmer;
		double dFPR;
		uint64_t nEntry;
		uint64_t tEntry;
	};

	BloomMap<T>(size_t filterSize, unsigned hashNum, unsigned kmerSize) :
			m_size(filterSize), m_hashNum(hashNum), m_dFPR(0), m_nEntry(0), m_tEntry(
					0), m_kmerSize(kmerSize) {
		m_array = new T[m_size]();
	}

	BloomMap<T>(const string &FilePath) {

		m_array = new T[m_size]();
		FILE *file = fopen(FilePath.c_str(), "rb");
		if (file == NULL) {
			cerr << "file \"" << FilePath << "\" could not be read." << endl;
			exit(1);
		}

		loadHeader(file);

		size_t countRead = fread(m_array, m_size, 1, file);
		if (countRead != 1 && fclose(file) != 0) {
			cerr << "file \"" << FilePath << "\" could not be read." << endl;
			exit(1);
		}
	}

	~BloomMap() {

		delete[] m_array;
	}

	T& operator[](size_t i) {
		assert(i < m_size);
		return m_array[i];
	}

	const T& operator[](size_t i) const {
		assert(i < m_size);
		return m_array[i];
	}

	void insert(std::vector<size_t> const &hashes, std::vector<T> &values) {
		assert(hashes.size() == m_hashNum);
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
#pragma omp atomic write
			m_array[pos] = values[i];
		}
	}

	std::vector<T> query(std::vector<size_t> const &hashes) {
		assert(hashes.size() == m_hashNum);
		std::vector<T> values(hashes.size());

		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
			values[i] = m_array[pos];

		}
		return values;
	}

// Calculating Pop Count
	size_t PopCnt() const {
		size_t i, popBF = 0;
		for (i = 0; i < m_size; ++i) {
			if (m_array[i] != 0) {
				++popBF;
			}
		}
		return popBF;
	}

// Calculating FPR based on the PopCount

	double getFPR() const {
		return pow(double(PopCnt()) / double(m_size), m_hashNum);
	}

// Calculate FPR based on hash functions, size and number of entries
// see http://en.wikipedia.org/wiki/Bloom_filter

	double calcFPR_numInserted(size_t numEntr) const {
		return pow(
				1.0
						- pow(1.0 - 1.0 / double(m_size),
								double(numEntr) * m_hashNum), double(m_hashNum));
	}

// Calculates the optimal FPR to use based on hash functions

	double calcFPR_hashNum(unsigned hashFunctNum) const {
		return pow(2, -hashFunctNum);
	}

	void loadHeader(FILE *file) {

		FileHeader header;

		if (fread(&header, sizeof(struct FileHeader), 1, file) == 1) {

			cerr << "Loading header..." << endl;

		} else {
			cerr << "Failed to header" << endl;
		}
		char magic[9];
		strncpy(magic, header.magic, 8);
		magic[8] = '\0';

		m_size = header.size;
		m_hashNum = header.nhash;
		m_kmerSize = header.kmer;
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
		header.nEntry = m_nEntry;
		header.tEntry = m_tEntry;
		out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));

	}

	void storeFilter(string const &filterFilePath) const {
		ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

		assert(myFile);

		writeHeader(myFile);

		//write out each block
		myFile.write((char*) m_array, m_size * sizeof(T));

		myFile.close();
		assert(myFile);
	}
private:

	size_t m_size;
	unsigned m_hashNum;
	T* m_array;
	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;
	unsigned m_kmerSize;
};

#endif /* BLOOMMAP_HPP_ */
