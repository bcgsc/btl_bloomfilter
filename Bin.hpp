#ifndef BIN_HPP // NOLINT(llvm-header-guard)
#define BIN_HPP

#include <cstddef>
#include <iostream>
#include <vector>

typedef uint64_t T;

// Forward declaraions.
class Bin;

class Bin
{
  public:
	Bin() = default;
	Bin(T& word, unsigned bitsPerBin, unsigned binNumber, T bitmask)
	  : m_word(word)
	  , m_bitsPerBin(bitsPerBin)
	  , m_binNumber(binNumber)
	  , m_bitmask(bitmask)
	{}
	T maxValue() const { return m_bitmask; };
	operator uint64_t() const { return (m_word >> (m_binNumber * m_bitsPerBin)) & m_bitmask; };
	Bin& operator=(const T newVal)
	{
		if (newVal > maxValue()) {
			std::cerr << "error: cannot set bin to values larger than " << maxValue() << "."
			          << std::endl;
		}
		T oldValMask = (maxValue() << (m_binNumber * m_bitsPerBin)) ^ ((T)0 - 1);
		m_word = (m_word & oldValMask) + (newVal << (m_binNumber * m_bitsPerBin));
		return *this; // returns this for chaining
	};
	Bin& operator=(const Bin& bin);

  private:
	/** Reference to the word containing this bin. */
	T& m_word;
	/** Number of bits in each counter. */
	unsigned m_bitsPerBin = 0;
	/** The nth bin this class is representing. */
	unsigned m_binNumber = 0;
	/** Masking bit used in bit operations. E.g. 0b 0000 0011 */
	T m_bitmask = 0;
};

Bin&
Bin::operator=(const Bin& bin)
{
	if (this == &bin) {
		return *this;
	}

	if (m_bitsPerBin != bin.m_bitsPerBin) {
		std::cerr << "Error: bitsPerBin do not match" << std::endl;
		exit(EXIT_FAILURE);
	}

	T newVal = (bin.m_word >> (bin.m_binNumber * bin.m_bitsPerBin)) & bin.m_bitmask;
	T oldValMask = (maxValue() << (m_binNumber * m_bitsPerBin)) ^ ((T)0 - 1);
	m_word = (m_word & oldValMask) + (newVal << (m_binNumber * m_bitsPerBin));

	// return the existing object so we can chain this operator
	return *this;
}

#endif // BIN_HPP
