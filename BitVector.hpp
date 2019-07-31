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
		return bytes * sizeof(T) / bitsPerCounter;
	}
	BitVector() = default;
	BitVector(size_t sz, unsigned bitsPerCounter)
	  : m_data(sz, 0)
	  , m_size(sz)
	  , m_bitsPerCounter(bitsPerCounter)
	  , m_numPartitions(sizeof(T) * 8 / bitsPerCounter)
	{
		unsigned* found = std::find(
		    std::begin(allowedBitsPerCounter), std::end(allowedBitsPerCounter), bitsPerCounter);
		if (found != std::end(allowedBitsPerCounter)) {
			m_maskingBits = ((T)1 << bitsPerCounter) - (T)1;
		} else {
			std::cerr << "ERROR: invalid bitsPerCounter value"
			          << "\n"
			          << "Accepted values are: 2, 4, 8, and 64" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	T operator[](size_t i) const
	{
		size_t pos = i / m_numPartitions;
		size_t sub_pos = i % m_numPartitions;
		return (m_data[pos] >> (sub_pos * m_bitsPerCounter)) & m_maskingBits;
	}
	bool atomicIncrement(size_t hash);
	// Conventional insertion function that calls atomicIncrement()
	void insert(size_t hash) { atomicIncrement(hash); };
	unsigned bitsPerCounter() const { return m_bitsPerCounter; };
	size_t size() const { return m_size; };
	T maxValue() const { return m_maskingBits; };
	size_t sizeInBytes() const { return m_size / sizeof(T) * m_bitsPerCounter; };
	const std::vector<T>& vector() { return m_data; };

  private:
	/** A vector of elements of type T. */
	std::vector<T> m_data;
	/** Size of vector (number of counters). */
	size_t m_size = 0;
	/** Masking bit used in bit operations. E.g. 0b 0000 0011 */
	T m_maskingBits = 0;
	/** Number of bits in each counter. */
	unsigned m_bitsPerCounter = 0;
	/** Number of partitions in each element of the vector. */
	unsigned m_numPartitions = 0;
};

inline bool
BitVector::atomicIncrement(size_t hash)
{
	size_t pos = hash / m_numPartitions;
	size_t sub_pos = hash % m_numPartitions;
	T oldWord = m_data[pos];
	T oldBits = (m_data[pos] >> (sub_pos * m_bitsPerCounter)) & m_maskingBits;
	if (oldBits == maxValue()) {
		return false;
	}
	T newWord = oldWord + ((T)1 << (sub_pos * m_bitsPerCounter));
	return __sync_bool_compare_and_swap(&m_data[pos], oldWord, newWord);
}

#endif // BitVector_HPP
