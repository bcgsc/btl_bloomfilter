/*
 * MIBFConstructSupport.hpp
 *
 * Purpose: To provide support for easier filter construction
 *
 * IDs are still managed by the class calling this object
 *
 * In order to use object one must have a iterator with same interface
 * as ntHashIterator
 *
 * Assumes saturation bit is being used
 * TODO: add functionality of strand bit
 *
 *  Created on: Mar 2028
 *      Author: Justin Chu
 */

#ifndef MIBFCONSTRUCTSUPPORT_HPP_
#define MIBFCONSTRUCTSUPPORT_HPP_

#include "MIBloomFilter.hpp"
#include <tuple>
#include <google/dense_hash_map>
#include <google/sparse_hash_map>
#include <google/dense_hash_set>
#include <sdsl/int_vector.hpp>

// T = T ID type, H = rolling hash itr
template<typename T, class H>
class MIBFConstructSupport {
public:
	/*
	 * numhashfunctions may also mean number of spaced seeds
	 */
	MIBFConstructSupport(size_t expectedEntries, unsigned k,
			unsigned numHashFunction, double occupancy,
			const vector<string> &spacedSeeds = vector<string>(0)) :
			m_isBVMade(false), m_isMIBFMade(false), m_expectedEntries(
					expectedEntries), m_k(k), m_h(numHashFunction), m_occupancy(
					occupancy), m_spacedSeeds(spacedSeeds), m_counts(
					vector<T>()) {
		m_filterSize = MIBloomFilter<T>::calcOptimalSize(m_expectedEntries,
				m_h, m_occupancy);
		m_bv = sdsl::bit_vector(m_filterSize);
	}

	~MIBFConstructSupport() {
		assert(m_isBVMade & m_isMIBFMade);
	}

	/*
	 * Returns count of collisions (counts unique k-mers)
	 */
	inline size_t insertBVColli(H &itr) {
		assert(!m_isBVMade);
		size_t count = 0;
		/* init rolling hash state and compute hash values for first k-mer */
		for (; itr != itr.end(); ++itr) {
			unsigned colliCount = 0;
			for (unsigned i = 0; i < m_h; ++i) {
				uint64_t pos = (*itr)[i] % m_bv.size();
				uint64_t *dataIndex = m_bv.data() + (pos >> 6);
				uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
				colliCount += __sync_fetch_and_or(dataIndex, bitMaskValue)
						>> (pos & 0x3F) & 1;
			}
			if (colliCount == m_h) {
				++count;
			}
		}
		return count;
	}

	void insertBV(H &itr) {
		assert(!m_isBVMade);
		/* init rolling hash state and compute hash values for first k-mer */
		for (; itr != itr.end(); ++itr) {
			for (unsigned i = 0; i < m_h; ++i) {
				uint64_t pos = (*itr)[i] % m_bv.size();
				uint64_t *dataIndex = m_bv.data() + (pos >> 6);
				uint64_t bitMaskValue = (uint64_t) 1 << (pos & 0x3F);
				(void) (__sync_fetch_and_or(dataIndex, bitMaskValue)
						>> (pos & 0x3F) & 1);
			}
		}
	}

	/*
	 * Generate empty miBF, can currently only be called once per object
	 */
	MIBloomFilter<T> *getEmptyMIBF() {
		assert(!m_isBVMade);
		m_isBVMade = true;
		MIBloomFilter<T> *miBF = new MIBloomFilter<T>(m_h, m_k, m_bv,
				m_spacedSeeds);
		m_counts = vector<T>(miBF->getPop(), 0);
		return miBF;
	}

	/*
	 * Uses single value Reservoir sampling
	 * pair<ID,ID> first ID stores the currentID, and ID stores the current observation count
	 * If the second ID exceeds max possible count, the ID is a critical ID
	 * Critical IDs are needed for partial hits and always replace existing IDs
	 *
	 * Once saturation is set, insertions are prevented
	 */
	//void insertMIBF(MIBloomFilter<T> &miBF, H &itr, T id) {
	void insertMIBF(MIBloomFilter<T> &miBF, H &itr, unsigned pos) {
		assert(m_isBVMade & !m_isMIBFMade);
		assert(m_isBVMade & !m_isMIBFMade);
		//get positions
		hashSet values;
		values.set_empty_key(miBF.size());
		
		unsigned base_pos = pos; // more accurate pos calculation
		while (itr != itr.end()) {
			for (unsigned i = 0; i < m_h; ++i) {
		//		values.insert((*itr)[i]);
		//	}
			//++itr;
		//}
		//for (hashSet::iterator itr = values.begin(); itr != values.end();
		//		itr++) {
				//uint64_t randomSeed = *itr ^ id;
				uint64_t randomSeed = (*itr)[i] ^ pos;
				bool strand = (itr).get_strand();  // get strand info. TODO
				//uint64_t rank = miBF.getRankPos(*itr);
				uint64_t rank = miBF.getRankPos((*itr)[i]);
				T count = __sync_add_and_fetch(&m_counts[rank], 1);
				T randomNum = std::hash<T> { }(randomSeed) % count;
				if (randomNum == count - 1) {
					//miBF.setData(rank, id);
					//miBF.setData(rank, pos);
					miBF.setData(rank, pos, strand); // send strand info TODO
				}	
			}
			++itr;
			//++pos; // not sure if pos is accurate
			pos = base_pos + itr.pos(); // more accurate pos calculation
		}
	}
		
	//void insertSaturation(MIBloomFilter<T> &miBF, H &itr, T id, unsigned start_pos) {
	void insertSaturation(MIBloomFilter<T> &miBF, H &itr, unsigned start_pos) {
		if (!m_isMIBFMade) {
			assert(m_isBVMade);
			m_isMIBFMade = true;
		}
		typedef google::dense_hash_set<uint64_t> SatSet;
		SatSet satVal;
		satVal.set_empty_key(miBF.size());
		//setSatIfMissing(miBF, id, itr);
		setSatIfMissing(miBF, itr, start_pos);
	}

	/*
	 * Returns number of bits in top level bit vector
	 */
	size_t getFilterSize() const {
		return m_filterSize;
	}

private:
	typedef google::dense_hash_set<uint64_t> hashSet;

	bool m_isBVMade;
	bool m_isMIBFMade;
	size_t m_expectedEntries;
	unsigned m_k, m_h;
	double m_occupancy;
	const vector<string> &m_spacedSeeds;
	sdsl::bit_vector m_bv;
	vector<T> m_counts;
	size_t m_filterSize;

	/*
	 * Attempts of mutate values to prevent saturation
	 * If unable to save values it will saturate values
	 * Small chance that mutation may erase entries
	 */
	//inline void setSatIfMissing(MIBloomFilter<T> &miBF, T id, H &itr, unsigned pos) {
	inline void setSatIfMissing(MIBloomFilter<T> &miBF, H &itr, unsigned pos) {
		unsigned base_pos = pos;
		while (itr != itr.end()) {
			//for each set of hash values, check for saturation
			vector<uint64_t> rankPos = miBF.getRankPos(*itr);
			vector<T> results = miBF.getData(rankPos);
			vector<T> replacementIDs(m_h);
			bool valueFound = false;
			vector<T> seenSet(m_h);
			for (unsigned i = 0; i < m_h; ++i) {
				T currentResult = results[i] & MIBloomFilter<T>::s_antiMask;
				//if (currentResult == id) {
				if (currentResult == pos){
					valueFound = true;
					break;
				}
				if (find(seenSet.begin(), seenSet.end(), currentResult)
						== seenSet.end()) {
					seenSet.push_back(currentResult);
				} else {
					replacementIDs.push_back(currentResult);
				}
			}
			if (!valueFound) {
				uint64_t replacementPos = m_counts.size();
				T minCount = numeric_limits<T>::min();
				for (unsigned i = 0; i < m_h; ++i) {
					T currentResult = results[i]
							& MIBloomFilter<T>::s_antiMask;
					if (find(replacementIDs.begin(), replacementIDs.end(),
							currentResult) != replacementIDs.end()) {
						if (minCount < m_counts[rankPos[i]]) {
							minCount = m_counts[rankPos[i]];
							replacementPos = rankPos[i];
						}
					}
				}
				//mutate if possible
				if (replacementPos != m_counts.size()) {
					//miBF.setData(replacementPos, id);
					miBF.setData(replacementPos, pos);
#pragma omp atomic update
					++m_counts[replacementPos];
				} else {
					miBF.saturate(*itr);
				}
			}
			++itr;
			pos = base_pos + itr.pos();

		}
	}
};
#endif /* MIBFCONSTRUCTSUPPORT_HPP_ */
