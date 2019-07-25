#ifndef Bin_HPP // NOLINT(llvm-header-guard)
#define Bin_HPP

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
	{
		m_value = (m_word >> (m_binNumber * m_bitsPerBin)) & m_bitmask;
	}
	T maxValue() const { return m_bitmask; };
	operator uint64_t() const { return m_value; };
	Bin& operator=(const T newVal)
	{
		if (newVal > maxValue()) {
			std::cerr << "error: cannot set bin to values larger than " << maxValue() << "."
			          << std::endl;
		}
		if (newVal > m_value) {
			m_word = m_word + ((newVal - m_value) << (m_binNumber * m_bitsPerBin));
		} else if (newVal < m_value) {
			m_word = m_word - ((m_value - newVal) << (m_binNumber * m_bitsPerBin));
		} else {
			return *this;
		}
		m_value = newVal;
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
	/** Value of this bin */
	T m_value = 0;
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

	if (bin.m_value > m_value) {
		m_word = m_word + ((bin.m_value - m_value) << (m_binNumber * m_bitsPerBin));
	} else if (bin.m_value < m_value) {
		m_word = m_word - ((m_value - bin.m_value) << (m_binNumber * m_bitsPerBin));
	} else {
		return *this;
	}
	m_value = bin.m_value;

	// return the existing object so we can chain this operator
	return *this;
}

#endif // Bit_HPP
