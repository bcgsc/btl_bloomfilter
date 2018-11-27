#ifndef COUNTBLOOM_H
#define COUNTBLOOM_H

#include <cstring>
#include <limits>
#include <cmath>
#include <map>

// Forward declaraions.
template <typename T>
class CountBloomFilter;

// Method declarations.
template <typename T>
std::ostream& operator<<(std::ostream&, const CountBloomFilter<T>&);

template <typename T>
class CountBloomFilter {
public:
	CountBloomFilter(): m_filter(NULL), m_size(0), m_sizeInBytes(0), m_hashNum(0),
			    m_kmerSize(0), m_nEntry(0), m_tEntry(0),
			    m_dFPR(0), m_countThreshold(0) { }
	CountBloomFilter(size_t sz, unsigned m_hashNum, unsigned m_kmerSize,
			 int m_countThreshold = 0)
		: m_hashNum(m_hashNum), m_kmerSize(m_kmerSize), m_nEntry(0), m_tEntry(0),
		  m_dFPR(0), m_countThreshold(m_countThreshold) {
		// Round up sz to a multiple of 64.
		m_size = ((sz - 1) | 63) + 1;
		m_filter  = new T[m_size];
		std::memset(m_filter, 0, sizeof(T) * m_size);
		m_sizeInBytes = m_size * sizeof(T);
	}
	CountBloomFilter(const string &path);
	~CountBloomFilter() {
		delete[] m_filter;
	}
	template <typename U> T      minCount(const U &hashes) const;
	template <typename U> size_t getMinCounter(const U &hashes) const;
	template <typename U> bool   contains(const U &hashes) const;
	template <typename U> void   insert(const U &hashes);
	template <typename U> bool   insertAndCheck(const U &hashes);
	template <typename U> void   incrementMin(const U &hashes);
	template <typename U> void   incrementAll(const U &hashes);
	unsigned getKmerSize(void) const;
	unsigned getHashNum(void) const;
	size_t   getFilterSize() const;
	size_t   popCount() const;
	double   FPR(void) const;

	// Serialization interface
	struct FileHeader {
		char magic[8];
		uint32_t hlen;
		uint64_t size;
		uint32_t nhash;
		uint32_t kmer;
		double   dFPR;
		uint64_t nEntry;
		uint64_t tEntry;
	};
	void readHeader(FILE *file);
	void readFilter(const string &path);
	void writeHeader(std::ostream& out) const;
	void writeFilter(string const &path) const;
	friend std::ostream& operator<< <> (std::ostream&, const CountBloomFilter&);
	
	// Debug utitlities.

	// A hash table to ever check filter's counts with accurate counts.
	map<string, uint64_t> htab;

private:
	// m_filter         : An array of elements of type T; the bit-array or filter.
	// m_size           : Size of bloom filter (size of m_filter array).
	// m_sizeInBytes    : Size of the bloom filter in bytes (m_size * sizeof(T)).
	// m_hashNum        : Number of hash functions.
	// m_kmerSize       : Size of a k-mer.
	// m_nEntry         : Number of items the bloom filter holds.
	// m_tEntry         : ?
	// m_dFPR           : Why d?
	// m_countThreshold : A count greater or equal to this threshold
	//                    establishes existence of an element in the filter.

	T        *m_filter;
	size_t   m_size;
	size_t   m_sizeInBytes;
	unsigned m_hashNum;
	unsigned m_kmerSize;
	size_t   m_nEntry;
	size_t   m_tEntry;
	double   m_dFPR;
	size_t   m_countThreshold;
};

// Method definitions.

// Bloom filter operations.

// Get the minimum count of the m_hashNum positions of m_filter.
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

// Get the index of the minimum counter from the m_hashNum positions of m_filter.
template <typename T>
template <typename U>
size_t CountBloomFilter<T>::getMinCounter(const U &hashes) const {
	size_t minPos = hashes[0] % m_size;
	T min = m_filter[minPos];
	for (size_t i = 1; i < m_hashNum; ++i) {
		size_t pos = hashes[i] % m_size;
		if (m_filter[pos] < min) {
			min = m_filter[pos];
			minPos = i;
		}
	}
	return minPos;
}

/*
Use of atomic increments in incrementMin() and incrementAll():

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

// Of the m_hashNum counters, increment only the minimum.
template <typename T>
template <typename U>
void CountBloomFilter<T>::incrementMin(const U &hashes) {
	static T currentVal, newVal;
	size_t pos = getMinCounter(hashes);
	do {
		currentVal = m_filter[pos];
		newVal     = currentVal + 1;
		if (newVal < currentVal)
			break;
	} while(!__sync_bool_compare_and_swap(&m_filter[pos], currentVal, newVal));
}

// Increment all the m_hashNum counters.
template <typename T>
template <typename U>
void CountBloomFilter<T>::incrementAll(const U &hashes) {
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

/*
Check if an element exists. If the minimum count at the m_hashNum positions of
m_filter is more than or equal to a predefined count threshold, then the element
is said to be present in the Bloom filter. count() therefore returns true when
this condition is satisfied, or else, false.
*/
template <typename T>
template <typename U>
bool CountBloomFilter<T>::contains(const U &hashes) const {
	return minCount(hashes) > m_countThreshold;
}

template <typename T>
template <typename U>
void CountBloomFilter<T>::insert(const U &hashes) {
	incrementMin(hashes);
}

template <typename T>
template <typename U>
bool CountBloomFilter<T>::insertAndCheck(const U &hashes) {
	bool found = contains(hashes);
	incrementMin(hashes);
	return found;
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

// Serialization interface.
template <typename T>
CountBloomFilter<T>::CountBloomFilter(const string &path) {
	if (m_filter != NULL)
		delete[] m_filter;
	m_filter = new T[m_size];
	readFilter(path);
}
template <typename T>
void CountBloomFilter<T>::readFilter(const string &path) {
	FILE *fp;
	if ((fp = fopen(path.c_str(), "rb")) == NULL) {
		cerr << "ERROR: Failed to open file: " << path << "\n";
		exit(1);
	}
	readHeader(fp);
	long int lCurPos = ftell(fp);
	fseek(fp, 0, 2);
	size_t arraySizeOnDisk = ftell(fp) - sizeof(struct FileHeader);
	fseek(fp, lCurPos, 0);
	if (arraySizeOnDisk != m_sizeInBytes) {
		cerr << "ERROR: File size of " << path << " ("
		     << arraySizeOnDisk << " bytes), "
		     << "does not match size read from its header ("
		     << m_sizeInBytes << " bytes).\n";
		exit(1);
	}

	size_t nread = fread(m_filter, arraySizeOnDisk, 1, fp);
	if (nread != 1 && fclose(fp) != 0) {
		cerr << "ERROR: The bit array could not be read from the file: "
		     << path << "\n";
		exit(1);
	}
}
template <typename T>
void CountBloomFilter<T>::readHeader(FILE *fp) {
	FileHeader header;
	char magic[9];
	if (fread(&header, sizeof(struct FileHeader), 1, fp) != 1) {
		cerr << "Failed to read header\n";
		exit(1);
	}
	memcpy(magic, header.magic, 8);
	magic[8] = '\0';
	m_size        = header.size;
	m_hashNum     = header.nhash;
	m_kmerSize    = header.kmer;
	m_sizeInBytes = m_size * sizeof(T);
}
template <typename T>
void CountBloomFilter<T>::writeFilter(string const &path) const {
	ofstream ofs(path.c_str(), ios::out | ios::binary);
	cerr << "Writing a " << m_sizeInBytes << " byte filter to a file on disk.\n";
	ofs << *this;
	ofs.close();
	assert(ofs);
}
template <typename T>
void CountBloomFilter<T>::writeHeader(std::ostream &out) const {
	FileHeader header;
	char magic[9];
	memcpy(header.magic, "BlOOMFXX", 8);
	memcpy(magic, header.magic, 8);
	magic[8] = '\0';
	header.hlen   = sizeof(struct FileHeader);
	header.size   = m_size;
	header.nhash  = m_hashNum;
	header.kmer   = m_kmerSize;
	header.dFPR   = m_dFPR;
	header.nEntry = m_nEntry;
	header.tEntry = m_tEntry;
	out.write(reinterpret_cast<char*>(&header), sizeof(struct FileHeader));
	assert(out);
}
// Serialize the bloom filter to a C++ stream */
template <typename T>
std::ostream& operator<<(std::ostream &os, const CountBloomFilter<T>& cbf) {
	assert(os);
	cbf.writeHeader(os);
	assert(os);
	os.write(reinterpret_cast<char*>(cbf.m_filter), cbf.m_sizeInBytes);
	assert(os);
	return os;
}

#endif // COUNTBLOOM_H
