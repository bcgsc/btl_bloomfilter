# bloomfilter

The BTL C/C++ Common Bloom filters for bioinformatics projects, as well as any APIs created for other programming languages.

# usage example (C++)

Fast Bloom filter loading using the rolling hash function.

```C++
#include "BloomFilter.hpp"
#include <vector>
#include <string>

using namespace std;

/** stores state between calls to rolling hash */
struct RollingHashState {
	/* seed hash value for current k-mer */
	uint64_t hash;
	/* seed hash value for reverse complement of current k-mer */
	uint64_t rcHash;
};

int main(int argc, char** argv)
{
	/* test sequence */
	const string seq = "TAGAATCACCCAAAGA";
	/* k-mer size */
	const unsigned k = 5;
	/* number of Bloom filter hash functions */
	const unsigned numHashes = 4;
	/* size of Bloom filter (in bits) */
	const unsigned size = 1000;
	/* hash values for current k-mer */
	vector<size_t> hashes;

	/* init Bloom filter */
	BloomFilter bloom(size, numHashes, k);

	/* init rolling hash state and compute hash values for first k-mer */
	RollingHashState state;
	string kmer0 = seq.substr(0,k);
	hashes = bloom.multiHash(kmer0.c_str(), state.hash, state.rcHash);

	/* load k-mers into Bloom filter using rolling hash */
	for (unsigned i = 1; i < seq.length() - k + 1; ++i) {
		/* "roll" hash values right to current k-mer */
		char charOut = seq[i - 1];
		char charIn = seq[i + k - 1];
		hashes = bloom.multiHash(state.hash, state.rcHash, charOut, charIn);
		/* insert current k-mer into Bloom filter */
		bloom.insert(hashes);
	}

	return 0;
}
```

# files

* `BloomFilter.hpp`: main Bloom filter class
* `rolling.h`: rolling hash function (required by `BloomFilter.hpp`)
* `Tests/Unit`: unit tests
* `Tests/AdHoc`: ad-hoc tests

# unit tests

The unit tests may be compiled and run with:

	$ ./autogen.sh
	$ ./configure
	$ make check

To see more detailed output for the individual tests, run the binaries in `Tests/Unit` from the command line. (The ad-hoc tests in `Tests/AdHoc` may also be run in this way.)

# acknowledgements

This projects uses:
* [CATCH](https://github.com/philsquared/Catch) unit test framework for C/C++
* rolling hash implementation by Hamid Mohamadi

# Bloom filter file format

The specification of the Bloom filter file format is as follows:

1. magic
  * Description: bf magic string
  * Type: char[8]
  * Value: BLOOMFXX
2. hlen
  * Description: length of the header text
  * Type: uint32_t
  * Value:
3. header
  * Description: Plain header text
  * Type: char[hlen]
  * Value:
    * size
      * Description: The size of Bloom filter
      * Type: char[8]
      * Value:
    * nhash
      * Description: number of hashes
      * Type: uint32_t
      * Value:
    * kmer [optional]
      * Description: k-mer size
      * Type: uint32_t
      * Value:
    * seed [optional]
      * Description: initial seeds for different hashes
      * Type: uint64_t[nhash]
      * Value: [0,1, ..., nhash-1]
    * dFPR [optional]
      * Description: desired false positve rate
      * Type: double
      * Value:
    * aFPR [optional]
      * Description: approximate false positive rate
      * Type: double
      * Value:
    * rFPR [optional]
      * Description: redundant false positive rate
      * Type: double
      * Value:
    * nEntry [optional]
      * Description: number of entries
      * Type: uint64_t
      * Value:
    * tEntry [optional]
      * Description: total number of entries
      * Type: uint64_t
      * Value:
4. filter
  * Description: Bloom filter content
  * Type: uchar[size]
  * Value:
