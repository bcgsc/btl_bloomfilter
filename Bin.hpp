#ifndef Bin_HPP // NOLINT(llvm-header-guard)
#define Bin_HPP

#include <algorithm>
#include <cstddef>
#include <iostream>
#include <vector>

using std::size_t;

unsigned allowedbitsPerBin[3] = { 2, 4, 8 };

typedef uint64_t T;

// Forward declaraions.
class Bin;

class Bin
{
  public:
	Bin() = default;
	Bin(T& word, unsigned bitsPerBin, unsigned binNumber)
	  : m_word(word)
	  , m_bitsPerBin(bitsPerBin)
	  , m_binNumber(binNumber)
	{
		unsigned* found =
		    std::find(std::begin(allowedbitsPerBin), std::end(allowedbitsPerBin), bitsPerBin);
		if (found != std::end(allowedbitsPerBin)) {
			m_maskingBits = ((T)1 << bitsPerBin) - (T)1;
		} else {
			std::cerr << "ERROR: invalid bitsPerBin value"
			          << "\n"
			          << "Accepted values are: 2, 4, and 8" << std::endl;
			exit(EXIT_FAILURE);
		}
	}
	T maxValue() const { return m_maskingBits; };
	operator uint64_t() const { return (m_word >> (m_binNumber * m_bitsPerBin)) & m_maskingBits; };
	T operator=(const T newVal)
	{
		m_word = m_word + (newVal << (m_binNumber * m_bitsPerBin));
		return newVal; // returns newVal for chaining
	};

  private:
	/** A vector of elements of type T. */
	T& m_word;
	/** Number of bits in each counter. */
	unsigned m_bitsPerBin = 0;
	/** The nth bin this class is representing. */
	unsigned m_binNumber = 0;
	/** Masking bit used in bit operations. E.g. 0b 0000 0011 */
	T m_maskingBits = 0;
};

#endif // Bit_HPP
