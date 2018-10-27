#ifndef COUNTBLOOM_H
#define COUNTBLOOM_H

#include <string>
#include <cstring>
#include <limits>
#include <vector>
#include <map>

#if _OPENMP
# include <omp.h>
#endif

using namespace std;

template <typename T>
class CountBloomFilter {
public:
	CountBloomFilter(): m_filter(NULL), m_size(0), m_sizeInBytes(0), m_hashNum(0),
			    m_kmerSize(0),  m_nEntry(0), m_countMax(0) { }
	CountBloomFilter(size_t sz, unsigned m_hashNum, unsigned m_kmerSize, int m_countMax)
		: m_hashNum(m_hashNum), m_kmerSize(m_kmerSize), m_nEntry(0),
		  m_countMax(m_countMax)  {
		// Round up to a multiple of 64 and initialize filter array.
		m_size = ((sz - 1) | 63) + 1;
		m_filter  = new T[m_size];
		std::memset(m_filter, 0, sizeof(T) * m_size);
		m_sizeInBytes = m_size * sizeof(T);

#ifdef _OPENMP
		m_locks.resize(10000);
		for (size_t i = 0; i < m_locks.size(); ++i)
			omp_init_nest_lock(&(m_locks.at(i)));
#endif
	}
	~CountBloomFilter() {
		delete[] m_filter;
#ifdef _OPENMP
		for (size_t i = 0; i < m_locks.size(); ++i)
			omp_destroy_nest_lock(&(m_locks.at(i)));
#endif
	}
	T minCount(const size_t hashes[]) const;
	template <typename U>
	bool contains(const U &hashes) const;
	template <typename U>
	void insert(const U &hashes);
	unsigned getKmerSize(void) const;
	unsigned getHashNum(void) const;
	size_t getFilterSize() const;
#ifdef _OPENMP
	void acquireLocks(const size_t hashes[]);
	void releaseLocks(const size_t hashes[]);
#endif
	size_t popCount() const;
	double FPR(void) const;

	// STUBS:
	CountBloomFilter(const string& path);
	void loadFilter(const string&);
	friend std::ostream& operator<<(std::ostream& out,
					const CountBloomFilter& bloom) {
		assert(out);
		// bloom.writeHeader(out);
		assert(out);
		out.write(reinterpret_cast<char*>(bloom.m_filter), bloom.m_sizeInBytes);
		assert(out);
		return out;
	}
	// STUBS.

	map<string, uint64_t> htab; // To ever compare bloom counts with true counts.

private:
#ifdef _OPENMP
	bool acquireLockNonBlocking(size_t hashVal);
	void acquireLock(size_t hashVal);
	void releaseLock(size_t hashVal);
	omp_nest_lock_t& getLock(size_t hashVal);

	mutable std::vector<omp_nest_lock_t> m_locks;
	T        *m_filter;     // A array of elements of type T; the bit-array or filter.
	size_t   m_size;        // Size of bloom filter (number of counters in m_filter).
	size_t   m_sizeInBytes; // Size of the bloom filter in bytes (m_size * sizeof(T)).
	unsigned m_hashNum;     // Number of hash functions.
	unsigned m_kmerSize;    // Length of a k-mer.
	size_t   m_nEntry;      // Number of items the bloom filter holds.
	T        m_countMax;    // Maximum value a T can hold; use to check for overflow.
#endif
};

// Method definitions.
	
template <typename T>
T CountBloomFilter<T>::minCount(const size_t hashes[]) const {
	T min = m_filter[hashes[0] % m_size];
	for (size_t i = 1; i < m_hashNum; ++i) {
		size_t bit = hashes[i] % m_size;
#ifdef _OPENMP
		acquireLock(hashes[i]);
		if (m_filter[bit] < min)
			min = m_filter[bit];
		releaseLock(hashes[i]);
#else
		if (m_filter[bit] < min)
			min = m_filter[bit];
#endif
	}
	return min;
}

template <typename T>
template <typename U>
bool CountBloomFilter<T>::contains(const U &hashes) const {
	return minCount(hashes) >= m_countMax;
}

template <typename T>
template <typename U>
void CountBloomFilter<T>::insert(const U &hashes) {
	for (size_t i = 0; i < m_hashNum; ++i) {
		size_t bit = hashes[i] % m_size;
#ifdef _OPENMP
		acquireLock(hashes[i]);
		if (m_filter[bit] < std::numeric_limits<T>::max())
			m_filter[bit]++;
		releaseLock(hashes[i]);
#else
		if (m_filter[bit] < std::numeric_limits<T>::max())
			m_filter[bit]++;
#endif
	}
}

template <typename T>
unsigned CountBloomFilter<T>::getKmerSize(void) const { return m_kmerSize; }

template <typename T>
unsigned CountBloomFilter<T>::getHashNum(void) const { return m_hashNum; }

template <typename T>
size_t CountBloomFilter<T>::getFilterSize() const { return m_size; }


// Utilities

#ifdef _OPENMP
// Acquire all the locks, or else, quit.
template <typename T>
void CountBloomFilter<T>::acquireLocks(const size_t hashes[]) {
	bool acquiredAllLocks = false;
	while (!acquiredAllLocks) {
		acquiredAllLocks = true;
		for (size_t i = 0; i < m_hashNum; ++i) {
			if (!acquireLockNonBlocking(hashes[i])) {
				acquiredAllLocks = false;
				for (size_t j = 0; j < i; ++j)
					releaseLock(hashes[j]);
				break;
			}
		}
	}
}

template <typename T>
void CountBloomFilter<T>::releaseLocks(const size_t hashes[]) {
	for (size_t i = 0; i < m_hashNum; ++i)
		releaseLock(hashes[i]);
}

template <typename T>
bool CountBloomFilter<T>::acquireLockNonBlocking(size_t hashVal) {
	return omp_test_nest_lock(&getLock(hashVal));
}

template <typename T>
void CountBloomFilter<T>::acquireLock(size_t hashVal) {
	omp_set_nest_lock(&getLock(hashVal));
}

template <typename T>
void CountBloomFilter<T>::releaseLock(size_t hashVal) {
	omp_unset_nest_lock(&getLock(hashVal));
}

template <typename T>
omp_nest_lock_t& CountBloomFilter<T>::getLock(size_t hashVal) {
	size_t lockIndex = (hashVal % m_size) % m_locks.size();
	return m_locks.at(lockIndex);
}

#endif // _OPENMP

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
	return pow((double)popCount() / (double)m_size, m_hashNum);
}

#endif // COUNTBLOOM_H
