#ifndef COUNTBLOOM_H
#define COUNTBLOOM_H

#include <cstring>
#include <limits>
#include <cmath>

template <typename T>
class CountBloomFilter {
public:
	CountBloomFilter(size_t sz, unsigned m_hashNum, unsigned m_kmerSize, int m_countThreshold)
		: m_hashNum(m_hashNum), m_kmerSize(m_kmerSize), m_nEntry(0),
		  m_countThreshold(m_countThreshold)  {
		// Round up sz to a multiple of 64.
		m_size = ((sz - 1) | 63) + 1;
		m_filter  = new T[m_size];
		std::memset(m_filter, 0, sizeof(T) * m_size);
		m_sizeInBytes = m_size * sizeof(T);
	}
	~CountBloomFilter() {
		delete[] m_filter;
	}
	template <typename U>
	T minCount(const U &hashes) const;
	template <typename U>
	bool contains(const U &hashes) const;
	template <typename U>
	void insert(const U &hashes);
	unsigned getKmerSize(void) const;
	unsigned getHashNum(void) const;
	size_t getFilterSize() const;
	size_t popCount() const;
	double FPR(void) const;

private:
	// m_filter         : A array of elements of type T; the bit-array or filter.
	// m_size           : Size of bloom filter (number of counters in m_filter).
	// m_sizeInBytes    : Size of the bloom filter in bytes (m_size * sizeof(T)).
	// m_hashNum        : Number of hash functions.
	// m_kmerSize       : Size of a k-mer.
	// m_nEntry         : Number of items the bloom filter holds.
	// m_countThreshold : A count greater or equal to this threshold
	//                    establishes existence of an element in the filter.

	T        *m_filter;        
	size_t   m_size;           
	size_t   m_sizeInBytes;    
	unsigned m_hashNum;        
	unsigned m_kmerSize;       
	size_t   m_nEntry;         
	T        m_countThreshold; 
};

// Method definitions.

/*  Return the minimum count of the m_hashNum positions of m_filter. */
template <typename T>
template <typename U>
T CountBloomFilter<T>::minCount(const U &hashes) const {
	T min = m_filter[hashes[0] % m_size];
	for (size_t i = 1; i < m_hashNum; ++i) {
		size_t pos = hashes[i] % m_size;
		if (m_filter[pos] < min)
			min = m_filter[pos];
	}
	return min;
}

/*
Check if an element exists. If the minimum count at the m_hashNum positions of
m_filter is more than or equal to a predefined count threshold, then the element
is said to be present in the Bloom filter. count() therefore returns true when
this condition is satisfied, or else, false.
*/
template <typename T>
template <typename U>
bool CountBloomFilter<T>::contains(const U &hashes) const {
	return minCount(hashes) >= m_countThreshold;
}

/*
Increment all the m_hashNum positions of the m_filter array by one.

A atomic compare-and-swap (CAS) operation increments m_filter[pos]. The CAS
operation takes a memory location and a value that the caller believes the
location currently stores. If the memory location still holds that value when
the atomic compare-and-swap executes, then a new value is stored and 'true'
returned; otherwise, memory is left unchanged and 'false' returned.

The value of m_filter[pos] may be changed by another thread between a read from
that memory location and a write to it. The CAS operation is called in a loop
until it succeeds, which ensures that a write does not happen if some other
thread has incremented the value between this thread's read and write.

Note that CAS operations suffer from the ABA problem.
*/
template <typename T>
template <typename U>
void CountBloomFilter<T>::insert(const U &hashes) {
	static T currentVal, newVal;
	for (size_t i = 0; i < m_hashNum; ++i) {
		size_t pos = hashes[i] % m_size;
		do {
			currentVal = m_filter[pos];
			newVal     = currentVal + 1;
			if (newVal < currentVal)
				break;
		} while(!__sync_bool_compare_and_swap(&m_filter[pos], currentVal, newVal));
	}
}

template <typename T>
unsigned CountBloomFilter<T>::getKmerSize(void) const { return m_kmerSize; }

template <typename T>
unsigned CountBloomFilter<T>::getHashNum(void) const { return m_hashNum; }

template <typename T>
size_t CountBloomFilter<T>::getFilterSize() const { return m_size; }

/* Count the number of non-zero counters. */
template <typename T>
size_t CountBloomFilter<T>::popCount() const {
	size_t count = 0;
	for (size_t i = 0; i < m_size; ++i) {
		if (m_filter[i] != 0)
			++count;
	}
	return count;
}

template <typename T>
double CountBloomFilter<T>::FPR(void) const {
	return std::pow((double)popCount() / (double)m_size, m_hashNum);
}

#endif // COUNTBLOOM_H
