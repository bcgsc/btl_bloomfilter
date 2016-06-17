/*
 * CountingBloomFilter.hpp
 *
 *  Created on: Dec 17, 2015
 *      Author: gjahesh
 */

#ifndef COUNTINGBLOOMFilter_HPP_
#define COUNTINGBLOOMFilter_HPP_

#include <vector>
#include <valarray>
#include <math.h>
#include <cassert>
#include <string>
#include <limits>
#include <utility>
#include "BloomMap.hpp"

#ifdef _OPENMP
#include <omp.h>
#endif

template<typename T>
class CountingBloomFilter {
public:

	/** Default constructor */
	CountingBloomFilter() : m_size(0), m_hashNum(0), m_bMinInc(true),
		m_bloomMap() {}

	/**
	 * Constructor
	 *
	 * @param filterSize number of element slots
	 * @param hashNum number of hash functions
	 * @param minInc If true, use "minimum increment" style incrementing,
	 * where only the minimum count(s) for an element are increased
	 * during an increment operation. If false, increment all counters
	 * for the element.
	 */
	CountingBloomFilter(size_t filterSize, unsigned hashNum,
		bool minInc = true) :
		m_size(filterSize), m_hashNum(hashNum), m_bMinInc(minInc),
		m_bloomMap(filterSize, hashNum, 0) {}

	/** Destructor */
	virtual ~CountingBloomFilter() {
	}

	/** Return the count of the single element*/
	T operator[](size_t i) const {
		return m_bloomMap[i];
	}

	/**
	 * Add the object to this counting multiset.
	 * If all values are the same update all
	 * If some values are larger only update smallest counts
	 *
	 * @param hashes array of hash values for the element.  The
	 * templated type ArrayT is used so that the caller may use either
	 * a std::vector or a plain old C array.
	 * @return std::pair<T,bool> where T is the count value
	 * before incrementing, and bool is true if the count was
	 * successfully incremented. The bool component will be false when
	 * the counter has already reached its max value (saturation).
	 */
	template <typename ArrayT>
	std::pair<T,bool> insert(const ArrayT& hashes) {

		m_bloomMap.getLocks(hashes);

		//check for which elements to update, basically holding the minimum
		//hash value i.e counter value.
		T minEle = (*this)[hashes];

		//saturate at max counter value (don't roll over to 0)
		if (minEle == std::numeric_limits<T>::max()) {
			m_bloomMap.releaseLocks(hashes);
			return std::make_pair(minEle, false);
		}

		for (unsigned int i = 0; i < m_hashNum; ++i) {
			size_t hashVal = hashes[i] % m_size;
			T val = m_bloomMap[hashVal];
			// if m_bMinInc is true, update only those elements that
			// have a minimum counter value
			if (!m_bMinInc || minEle == val) {
				insert(hashVal);
			}
		}

		m_bloomMap.releaseLocks(hashes);
		return std::make_pair(minEle, true);
	}

	/** Add the object with the specified index (debug). */
	void insert(size_t index) {
#pragma omp atomic
		++m_bloomMap[index];
	}

	/**
	 * Return the count of an element based on a Minimum Selection
	 *
	 * @param hashes array of hash values for the element.  The
	 * templated type ArrayT is used so that the caller may use either
	 * a std::vector or a plain old C array.
	 */
	template <typename ArrayT>
	T operator[](const ArrayT& hashes) const {

		T currentMin = m_bloomMap[hashes[0] % m_size];
		for (unsigned int i = 1; i < m_hashNum; ++i) {
			T min = m_bloomMap[hashes[i] % m_size];
			if (min < currentMin) {
				currentMin = min;
			}
			if (0 == currentMin) {
				return (0);
			}
		}
		return currentMin;
	}

	/** Return the count of the single element*/

	T query(size_t i) const {
		return m_bloomMap[i];
	}

	/** Get size of Bloom filter (number of elements of type T) */
	size_t getFilterSize() const
	{
		return m_size;
	}

	/** Get Bloom filter false positive rate */
	double getFPR() const
	{
		return m_bloomMap.getFPR();
	}

private:

	size_t m_size;
	unsigned m_hashNum;
	bool m_bMinInc;
	BloomMap<T> m_bloomMap;
};

#endif /* COUNTINGBLOOM_HPP_ */
