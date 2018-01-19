/*
 * MIBloomFilter.hpp
 *
 * Agnostic of hash function used -> cannot call contains without an array of hash values
 *
 *  Created on: Jan 14, 2016
 *      Author: cjustin
 */

#ifndef MIBLOOMFILTER_HPP_
#define MIBLOOMFILTER_HPP_

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
#include <limits>
#include <sdsl/bit_vector_il.hpp>
#include <sdsl/rank_support.hpp>
#include <omp.h>

using namespace std;
template<typename T>
class MIBloomFilter {
public:
	static const unsigned BLOCKSIZE = 512;
	//static methods
	/*
	 * Parses spaced seed string (string consisting of 1s and 0s) to vector
	 */
	static vector<vector<unsigned> > parseSeedString(
			const vector<string> &spacedSeeds) {
		vector<vector<unsigned> > seeds(spacedSeeds.size(), vector<unsigned>());
		for (unsigned i = 0; i < spacedSeeds.size(); ++i) {
			const string ss = spacedSeeds.at(i);
			for (unsigned j = 0; j < ss.size(); ++j) {
				if (ss.at(j) == '0') {
					seeds[i].push_back(j);
				}
			}
		}
		return seeds;
	}

	/*
	 * Returns an a filter size large enough to maintain an occupancy specified
	 */
	static size_t calcOptimalSize(size_t entries, unsigned hashNum,
			double occupancy) {
		size_t non64ApproxVal = size_t(
				-double(entries) * double(hashNum) / log(1.0 - occupancy));
		return non64ApproxVal + (64 - non64ApproxVal % 64);
	}

	/*
	 * Inserts a set of hash values into an sdsl bitvector and returns the number of collisions
	 * Thread safe on the bv, though return values will not be the same run to run
	 */
	static unsigned insert(sdsl::bit_vector &bv, uint64_t * hashValues,
			unsigned hashNum) {
		unsigned colliCount = 0;
		for (size_t i = 0; i < hashNum; ++i) {
			size_t pos = hashValues[i] % bv.size();
			uint64_t *dataIndex = bv.data() + (pos >> 6);
			uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
			colliCount += __sync_fetch_and_or(dataIndex, bitMaskValue)
					>> (pos & 0x3F) & 1;
		}
		return colliCount;
	}

	//TODO: include allowed miss in header
#pragma pack(1) //to maintain consistent values across platforms
	struct FileHeader {
		char magic[8];
		uint32_t hlen;	//header length (including spaced seeds)
		uint64_t size;
		uint64_t nEntry;
		uint64_t tEntry;
		double dFPR;
		uint32_t nhash;
		uint32_t kmer;
//		uint8_t allowedMiss;
	};

	/*
	 * Constructor using a prebuilt bitvector
	 */
	MIBloomFilter<T>(size_t expectedElemNum, double fpr, unsigned hashNum,
			unsigned kmerSize, sdsl::bit_vector &bv, size_t unique,
			const vector<string> seeds = vector<string>(0)) :
			m_dSize(0), m_dFPR(fpr), m_nEntry(unique), m_tEntry(
					expectedElemNum), m_hashNum(hashNum), m_kmerSize(kmerSize), m_sseeds(
					seeds) {
		cerr << "Converting bit vector to rank interleaved form" << endl;
		double start_time = omp_get_wtime();
		m_bv = sdsl::bit_vector_il < BLOCKSIZE > (bv);
		bv = sdsl::bit_vector();
		double time = omp_get_wtime() - start_time;
		cerr << "Converted bit vector to rank interleaved form " << time << "s"
				<< endl;
		if (!seeds.empty()) {
			m_ssVal = parseSeedString(m_sseeds);
			assert(m_sseeds[0].size() == kmerSize);
			for (vector<string>::const_iterator itr = m_sseeds.begin();
					itr != m_sseeds.end(); ++itr) {
				//check if spaced seeds are all the same length
				assert(m_kmerSize == itr->size());
			}
		}
		m_rankSupport = sdsl::rank_support_il < 1 > (&m_bv);
		m_dSize = getPop();
		m_data = new T[m_dSize]();
	}

	MIBloomFilter<T>(const string &filterFilePath) {
#pragma omp parallel for
		for (unsigned i = 0; i < 2; ++i) {
			if (i == 0) {
				FILE *file = fopen(filterFilePath.c_str(), "rb");
				if (file == NULL) {
#pragma omp critical(stderr)
					cerr << "file \"" << filterFilePath
							<< "\" could not be read." << endl;
					exit(1);
				}

				FileHeader header;
				if (fread(&header, sizeof(struct FileHeader), 1, file) == 1) {
#pragma omp critical(stderr)
					cerr << "Loading header..." << endl;
				} else {
#pragma omp critical(stderr)
					cerr << "Failed to Load header" << endl;
					exit(1);
				}
				char magic[9];
				strncpy(magic, header.magic, 8);
				magic[8] = '\0';
#pragma omp critical(stderr)
				cerr << "Loaded header... magic: " << magic << " hlen: "
						<< header.hlen << " size: " << header.size << " nhash: "
						<< header.nhash << " kmer: " << header.kmer << " dFPR: "
						<< header.dFPR << " nEntry: " << header.nEntry
						<< " tEntry: " << header.tEntry << endl;

				m_dFPR = header.dFPR;
				m_nEntry = header.nEntry;
				m_hashNum = header.nhash;
				m_tEntry = header.tEntry;
				m_kmerSize = header.kmer;
				m_dSize = header.size;
				m_data = new T[m_dSize]();

				if (header.hlen > sizeof(struct FileHeader)) {
					//load seeds
					for (unsigned i = 0; i < header.nhash; ++i) {
						char temp[header.kmer];

						if (fread(temp, header.kmer, 1, file) != 1) {
							cerr << "Failed to load spaced seed string" << endl;
							exit(1);
						} else {
							cerr << "Spaced Seed " << i << ": "
									<< string(temp, header.kmer) << endl;
						}
						m_sseeds.push_back(string(temp, header.kmer));
					}

					m_ssVal = parseSeedString(m_sseeds);
					assert(m_sseeds[0].size() == m_kmerSize);
					for (vector<string>::const_iterator itr = m_sseeds.begin();
							itr != m_sseeds.end(); ++itr) {
						//check if spaced seeds are all the same length
						assert(m_kmerSize == itr->size());
					}
				}

#pragma omp critical(stderr)
				cerr << "Loading data vector" << endl;

				long int lCurPos = ftell(file);
				fseek(file, 0, 2);
				size_t fileSize = ftell(file) - header.hlen;
				fseek(file, lCurPos, 0);
				if (fileSize != m_dSize * sizeof(T)) {
					cerr << "Error: " << filterFilePath
							<< " does not match size given by its header. Size: "
							<< fileSize << " vs " << m_dSize * sizeof(T)
							<< " bytes." << endl;
					exit(1);
				}

				size_t countRead = fread(m_data, fileSize, 1, file);
				if (countRead != 1 && fclose(file) != 0) {
					cerr << "file \"" << filterFilePath
							<< "\" could not be read." << endl;
					exit(1);
				}
			} else {
				string bvFilename = filterFilePath + ".sdsl";
#pragma omp critical(stderr)
				cerr << "Loading sdsl interleaved bit vector from: "
						<< bvFilename << endl;
				load_from_file(m_bv, bvFilename);
				m_rankSupport = sdsl::rank_support_il < 1 > (&m_bv);
			}
		}

		cerr << "Bit Vector Size: " << m_bv.size() << endl;
		cerr << "Popcount: " << getPop() << endl;
	}

	/*
	 * Stores the filter as a binary file to the path specified
	 * Stores uncompressed because the random data tends to
	 * compress poorly anyway
	 */
	inline void store(string const &filterFilePath) const {

#pragma omp parallel for
		for (unsigned i = 0; i < 2; ++i) {
			if (i == 0) {
				ofstream myFile(filterFilePath.c_str(), ios::out | ios::binary);

				assert(myFile);
				writeHeader(myFile);

				cerr << "Storing filter. Filter is " << m_dSize * sizeof(T)
						<< " bytes." << endl;

				//write out each block
				myFile.write(reinterpret_cast<char*>(m_data),
						m_dSize * sizeof(T));

				myFile.close();
				assert(myFile);

				FILE *file = fopen(filterFilePath.c_str(), "rb");
				if (file == NULL) {
					cerr << "file \"" << filterFilePath
							<< "\" could not be read." << endl;
					exit(1);
				}
			} else {
				string bvFilename = filterFilePath + ".sdsl";

				cerr << "Storing sdsl interleaved bit vector to: " << bvFilename
						<< endl;
				store_to_file(m_bv, bvFilename);
				cerr << "Number of bit vector buckets is " << m_bv.size()
						<< endl;
				cerr << "Uncompressed bit vector size is "
						<< (m_bv.size() + m_bv.size() * 64 / BLOCKSIZE) / 8
						<< " bytes" << endl;
			}
		}
	}

	/*
	 * saturated should be set to true to start with, if already saturated it should return as false
	 */
	inline bool insert(const size_t *hashes, T value, unsigned max,
			bool &saturated) {
		bool someValueSet = false;
		unsigned count = 0;
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = m_rankSupport(hashes[i] % m_bv.size());
			//check for saturation
			T oldVal = setVal(&m_data[pos], value);
			if (oldVal > mask) {
				oldVal = oldVal & antiMask;
			} else {
				saturated = false;
			}
			if (oldVal == 0 || oldVal == value) {
				someValueSet = true;
				++count;
			}
			if (count >= max) {
				return someValueSet;
			}
		}
		//TODO call manually?
		if (!someValueSet && !saturated) {
			saturate(hashes);
		}
		return someValueSet;
	}

	inline void saturate(const size_t *hashes) {
		for (size_t i = 0; i < m_hashNum; ++i) {
			size_t pos = m_rankSupport(hashes[i] % m_bv.size());
			__sync_or_and_fetch(&m_data[pos], mask);
		}
	}

	inline vector<T> at(const size_t *hashes, bool &saturated,
			unsigned maxMiss = 0) {
		vector<T> results(m_hashNum);
		unsigned misses = 0;
		for (unsigned i = 0; i < m_hashNum; ++i) {
			size_t pos = hashes[i] % m_bv.size();
			if (m_bv[pos] == 0) {
				++misses;
				saturated = false;
				if (misses > maxMiss) {
					return vector<T>();
				}
			} else {
				size_t rankPos = m_rankSupport(pos);
				T tempResult = m_data[rankPos];
				if (tempResult > mask) {
					results[i] = m_data[rankPos] & antiMask;
				} else {
					results[i] = m_data[rankPos];
					saturated = false;
				}
			}
		}
		return results;
	}

	inline const vector<vector<unsigned> > &getSeedValues() const {
		return m_ssVal;
	}

	inline unsigned getKmerSize() {
		return m_kmerSize;
	}

	inline unsigned getHashNum() {
		return m_hashNum;
	}

	/*
	 * computes id frequency based on datavector
	 */
	inline void getIDCounts(vector<size_t> &counts) {
		for (size_t i = 0; i < m_dSize; ++i) {
			++counts[m_data[i] & antiMask];
		}
	}

	/*
	 * Return FPR based on popcount
	 */
	inline double getFPR() const {
		return pow(double(getPop()) / double(m_bv.size()), double(m_hashNum));
	}

	/*
	 * Return FPR based on popcount and minimum number of matches for a hit
	 */
	inline double getFPR(unsigned allowedMiss) const {
		assert(allowedMiss < m_hashNum);
		double cumulativeProb = 0;
		double popCount = getPop();
		double p = popCount / double(m_bv.size());
		for (unsigned i = m_hashNum - allowedMiss; i <= m_hashNum; ++i) {
			cumulativeProb += double(nChoosek(m_hashNum, i)) * pow(p, i)
					* pow(1.0 - p, (m_hashNum - i));
		}
		return (cumulativeProb);
	}

	/*
	 * Return FPR based on number of inserted elements
	 */
	inline double getFPR_numEle() const {
		assert(m_nEntry > 0);
		return calcFPR_numInserted(m_nEntry);
	}

	inline size_t getPop() const {
		size_t index = m_bv.size() - 1;
		while (m_bv[index] == 0) {
			--index;
		}
		return m_rankSupport(index - 1) + 1;
	}

	inline size_t getPopNonZero() const {
		size_t count = 0;
		for (size_t i = 0; i < m_dSize; ++i) {
			if (m_data[i] != 0) {
				++count;
			}
		}
		return count;
	}

	inline size_t getPopSaturated() const {
		size_t count = 0;
		for (size_t i = 0; i < m_dSize; ++i) {
			if (m_data[i] >= mask) {
				++count;
			}
		}
		return count;
	}

	inline size_t getUniqueEntries() const {
		return m_nEntry;
	}

	inline size_t size() {
		return m_bv.size();
	}

	~MIBloomFilter() {
		delete[] m_data;
	}

private:
	const T mask = 1 << (sizeof(T) * 8 - 1);
	const T antiMask = ~mask;

	/*
	 * Helper function for header storage
	 */
	inline void writeHeader(ofstream &out) const {
		FileHeader header;
		char magic[9];
		strncpy(magic, MAGICSTR, 8);
		magic[8] = '\0';
		strncpy(header.magic, magic, 8);

		header.hlen = sizeof(struct FileHeader) + m_kmerSize * m_sseeds.size();
		header.kmer = m_kmerSize;
		header.size = m_dSize;
		header.nhash = m_hashNum;
		header.dFPR = m_dFPR;
		header.nEntry = m_nEntry;
		header.tEntry = m_tEntry;

		cerr << "Writing header... magic: " << magic << " hlen: " << header.hlen
				<< " nhash: " << header.nhash << " size: " << header.size
				<< " dFPR: " << header.dFPR << " nEntry: " << header.nEntry
				<< " tEntry: " << header.tEntry << endl;

		out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));

		for (vector<string>::const_iterator itr = m_sseeds.begin();
				itr != m_sseeds.end(); ++itr) {
			out.write(itr->c_str(), m_kmerSize);
		}
	}

	/*
	 * Calculates the optimal number of hash function to use
	 * Calculation assumes optimal ratio of bytes per entry given a fpr
	 */
	inline static unsigned calcOptiHashNum(double fpr) {
		return unsigned(-log(fpr) / log(2));
	}

	/*
	 * Calculate FPR based on hash functions, size and number of entries
	 * see http://en.wikipedia.org/wiki/Bloom_filter
	 */
	inline double calcFPR_numInserted(size_t numEntr) const {
		return pow(
				1.0
						- pow(1.0 - 1.0 / double(m_bv.size()),
								double(numEntr) * double(m_hashNum)),
				double(m_hashNum));
	}

	/*
	 * Calculates the optimal FPR to use based on hash functions
	 */
	inline double calcFPR_hashNum(int hashFunctNum) const {
		return pow(2.0, -hashFunctNum);
	}

	/*
	 * Returns old value that was inside
	 */
	inline T setVal(T *val, T newVal) {
		T oldValue;
		do {
			oldValue = *val;
			if (oldValue != 0)
				break;
		} while (!__sync_bool_compare_and_swap(val, oldValue, newVal));
		return oldValue;
	}

	inline unsigned nChoosek(unsigned n, unsigned k) const {
		if (k > n)
			return 0;
		if (k * 2 > n)
			k = n - k;
		if (k == 0)
			return 1;

		int result = n;
		for (unsigned i = 2; i <= k; ++i) {
			result *= (n - i + 1);
			result /= i;
		}
		return result;
	}

	//size of bitvector
	size_t m_dSize;

	sdsl::bit_vector_il<BLOCKSIZE> m_bv;
	T* m_data;
	sdsl::rank_support_il<1> m_rankSupport;

	double m_dFPR;
	uint64_t m_nEntry;
	uint64_t m_tEntry;
	unsigned m_hashNum;
	unsigned m_kmerSize;

	typedef vector<vector<unsigned> > SeedVal;
	vector<string> m_sseeds;
	SeedVal m_ssVal;
	const char* MAGICSTR = "MIBLOOMF";
};

#endif /* MIBLOOMFILTER_HPP_ */
