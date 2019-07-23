#ifndef BitVector_HPP // NOLINT(llvm-header-guard)
#define BitVector_HPP

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

using std::size_t;

unsigned allowedBitsPerCounter[4] = { 2, 4, 8, 64 };

typedef uint64_t T;

// Forward declaraions.
class BitVector;

class BitVector
{
  public:
	static inline size_t bytesToElements(size_t bytes, unsigned bitsPerCounter)
	{
		return bytes / sizeof(T) * 8 / bitsPerCounter;
	}
	BitVector() = default;
	BitVector(size_t sz, unsigned bitsPerCounter)
	  : m_vector((sz / sizeof(T)), 0)
	  , m_size(sz / sizeof(T) * 8 / bitsPerCounter)
	  , m_sizeInBytes(sz)
	  , m_bitsPerCounter(bitsPerCounter)
	  , m_numPartitions(sizeof(T) * 8 / bitsPerCounter)
	{
		unsigned* found = std::find(
		    std::begin(allowedBitsPerCounter), std::end(allowedBitsPerCounter), bitsPerCounter);
		if (found != std::end(allowedBitsPerCounter)) {
			m_maskingBits = (1 << bitsPerCounter) - 1;
		} else {
			std::cerr << "ERROR: invalid bitsPerCounter value"
			          << "\n"
			          << "Accepted values are: 2, 4, 8, and 64" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	T operator[](size_t i)
	{
		size_t pos = i / m_numPartitions;
		size_t sub_pos = i % m_numPartitions;
		return (m_vector[pos] >> (sub_pos * m_bitsPerCounter)) & m_maskingBits;
	}
	bool atomicIncrement(size_t hash);
	// Conventional insertion function that calls atomicIncrement()
	void insert(size_t hash) { atomicIncrement(hash); };
	unsigned bitsPerCounter() const { return m_bitsPerCounter; };
	size_t size() const { return m_size; };
	size_t maxValue() const { return m_maskingBits; };
	size_t sizeInBytes() const { return m_sizeInBytes; };
	std::vector<T> vector() { return m_vector; };

  private:
	// m_vector             : A vector of elements of type T.
	// m_size               : Size of vector (number of counters).
	// m_sizeInBytes        : Size of the vector in bytes.
	// m_maskingBits        : Masking bit used in bit operations. E.g. 0b 0000 0011
	// m_bitsPerCounter     : Number of bits in each counter
	// m_numPartitions      : Number of partitions in each element of the vector
	// m_increment          : Increment Value

	std::vector<T> m_vector;
	size_t m_size = 0;
	size_t m_sizeInBytes = 0;
	size_t m_maskingBits = 0;
	unsigned m_bitsPerCounter = 0;
	unsigned m_numPartitions = 0;
	size_t m_incrementUnit = 1;
};

inline bool
BitVector::atomicIncrement(size_t hash)
{
	size_t pos = hash / m_numPartitions;
	size_t sub_pos = hash % m_numPartitions;
	size_t oldByte = m_vector[pos];
	size_t oldBits = (m_vector[pos] >> (sub_pos * m_bitsPerCounter)) & m_maskingBits;
	if (oldBits == maxValue()) {
		return false;
	}
	size_t newByte = oldByte + (m_incrementUnit << (sub_pos * m_bitsPerCounter));
	return __sync_bool_compare_and_swap(&m_vector[pos], oldByte, newByte);
}

#endif // BitVector_HPP
