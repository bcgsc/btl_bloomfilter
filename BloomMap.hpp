/*
 * BloomMap.hpp
 *
 *  Created on: Dec 17, 2015
 *      Author: gjahesh
 */

#ifndef BLOOMMAP_HPP_
#define BLOOMMAP_HPP_

#include <vector>
#include <fstream>
#include <iostream>
#include <stdio.h>

#ifdef _OPENMP
#include <omp.h>
#endif

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

	/** Default constructor */
	BloomMap<T>() : m_size(0), m_hashNum(0), m_array(NULL), m_dFPR(0),
		m_nEntry(0), m_tEntry(0), m_kmerSize(0) {}

	BloomMap<T>(size_t filterSize, unsigned hashNum, unsigned kmerSize) :
		m_size(filterSize), m_hashNum(hashNum), m_dFPR(0), m_nEntry(0),
		m_tEntry(0), m_kmerSize(kmerSize)
	{
		m_array = new T[m_size]();

#ifdef _OPENMP
		/* locks for ensuring thread-safety */
		m_locks.resize(10000);
		for (size_t i = 0; i < m_locks.size(); i++)
			omp_init_nest_lock(&(m_locks.at(i)));
#endif
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

		if (m_array != NULL)
			delete[] m_array;
#ifdef _OPENMP
		for (size_t i = 0; i < m_locks.size(); i++)
			omp_destroy_nest_lock(&(m_locks.at(i)));
#endif
	}

	T& operator[](size_t i) {
		assert(i < m_size);
		assert(m_array != NULL);
		return m_array[i];
	}

	const T& operator[](size_t i) const {
		assert(i < m_size);
		assert(m_array != NULL);
		return m_array[i];
	}

	void insert(std::vector<size_t> const &hashes, std::vector<T> &values) {
		assert(hashes.size() == m_hashNum);
		assert(m_array != NULL);
		//iterates through hashed values adding it to the filter
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
#ifdef _OPENMP
			getLock(pos);
			m_array[pos] = values[i];
			releaseLock(pos);
#else
			m_array[pos] = values[i];
#endif
		}
	}

	std::vector<T> query(std::vector<size_t> const &hashes) {
		assert(hashes.size() == m_hashNum);
		assert(m_array != NULL);
		std::vector<T> values(hashes.size());

		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes.at(i) % m_size;
			assert(pos < m_size);
#ifdef _OPENMP
			getLock(pos);
			values[i] = m_array[pos];
			releaseLock(pos);
#else
			values[i] = m_array[pos];
#endif
		}
		return values;
	}

// Calculating Pop Count
	size_t PopCnt() const {
		assert(m_array != NULL);
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
		assert(m_array != NULL);
		ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

		assert(myFile);

		writeHeader(myFile);

		//write out each block
		myFile.write((char*) m_array, m_size * sizeof(T));

		myFile.close();
		assert(myFile);
	}

#ifdef _OPENMP
	/**
	 * Acquire all locks for the Bloom filter positions associated
	 * with an element. In order to avoid deadlocks, the locks must
	 * either be acquired in a well-defined order or a thread must try
	 * to acquire all of its locks at the same time and "back off" when
	 * this is not possible. I have opted for the latter approach, as sorting
	 * the hash values to produce a predictable lock ordering would
	 * be too expensive.
	 *
	 * @param hashes array of hash values for the element.  The
	 * templated type ArrayT is used so that the caller may use either
	 * a std::vector or a plain old C array.
	 */
	template <typename ArrayT>
	void getLocks(const ArrayT& hashes)
	{
		bool acquiredAllLocks = false;
		while (!acquiredAllLocks) {
			acquiredAllLocks = true;
			for (unsigned i = 0; i < m_hashNum; ++i) {
				size_t iLockIndex = hashes[i] % m_size % m_locks.size();
				assert(iLockIndex < m_locks.size());
				/**
				 * tricky: a nested lock must be used here because two
				 * distinct hash values for the same element may map to
				 * the same lock index
				 */
				if (!omp_test_nest_lock(&(m_locks.at(iLockIndex)))) {
					acquiredAllLocks = false;
					for (unsigned j = 0; j < i; ++j) {
						size_t jLockIndex = hashes[j] % m_size % m_locks.size();
						assert(jLockIndex < m_locks.size());
						omp_unset_nest_lock(&(m_locks.at(jLockIndex)));
					}
					break;
				}
			}
		}
	}

	/**
	 * Release all locks for the Bloom filter positions associated
	 * with an element.
	 *
	 * @param hashes array of hash values for the element.  The
	 * templated type ArrayT is used so that the caller may use either
	 * a std::vector or a plain old C array.
	 */
	template <typename ArrayT>
	void releaseLocks(const ArrayT& hashes)
	{
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t lockIndex = hashes[i] % m_size % m_locks.size();
			assert(lockIndex < m_locks.size());
			omp_unset_nest_lock(&(m_locks.at(lockIndex)));
		}
	}
#endif

private:

#ifdef _OPENMP

	void getLock(size_t pos)
	{
		assert(pos < m_size);
		size_t lockIndex = pos % m_locks.size();
		omp_set_nest_lock(&(m_locks.at(lockIndex)));
	}

	void releaseLock(size_t pos)
	{
		assert(pos < m_size);
		size_t lockIndex = pos % m_locks.size();
		omp_unset_nest_lock(&(m_locks.at(lockIndex)));
	}
#endif

	size_t m_size;
	unsigned m_hashNum;
	T* m_array;
	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;
	unsigned m_kmerSize;

#ifdef _OPENMP
	std::vector<omp_nest_lock_t> m_locks;
#endif
};

#endif /* BLOOMMAP_HPP_ */
