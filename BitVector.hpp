#ifndef BitVector_HPP // NOLINT(llvm-header-guard)
#define BitVector_HPP

#include <cstddef>
#include <vector>
#include <iostream>

using std::size_t;

// Forward declaraions.
template<typename T>
class BitVector;

template<typename T>
class BitVector
{
  public:
	BitVector() = default;
	BitVector(size_t sz, unsigned bitsPerCounter)
	  : m_vector((sz / sizeof(T)), 0)
	  , m_size(sz / sizeof(T) * 8 / bitsPerCounter)
	  , m_sizeInBytes(sz)
	  , m_bitsPerCounter(bitsPerCounter)
	  , m_numPartitions(sizeof(T) * 8 / bitsPerCounter)
	{
		switch (bitsPerCounter) {
		case 2:
			m_maskingBits = 3; //   equivalent to 0b 0000 0011 if T = uint8_t
			break;
		case 4:
			m_maskingBits = 15; //  equivalent to 0b 0000 1111 if T = uint8_t
			break;
		case 8:
			m_maskingBits = 255; // equivalent to 0b 1111 1111 if T = uint8_t
			break;
		default:
			std::cerr << "ERROR: invalid bitsPerCounter value"
			          << "\n"
			          << "Accepted values are: 2, 4 and 8" << std::endl;
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

template<typename T>
inline bool
BitVector<T>::atomicIncrement(size_t hash)
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
