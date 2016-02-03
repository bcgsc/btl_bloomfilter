# bloomfilter

The BTL C/C++ Common Bloom filters for bioinformatics projects, as well as any APIs created for other programming languages.

# usage example (C++)

Fast Bloom filter loading using the rolling hash function.

```C++
#include "BloomFilter.hpp"
#include <vector>
#include <string>
#include "RollingHashIterator.h"

using namespace std;

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
	RollingHashIterator itr(seq, numHashes, k);
	while (itr != itr.end()) {
		BloomFilterFilter.insert(*itr);
		itr++;
	}

	return 0;
}
```

# files

* `BloomFilter.hpp`: main Bloom filter class
* `RollingHashIterator.h`: Enable rolling hashing on a string 
* `RollingHash.h`: rolling hash interface (required by `RollingHashIterator.h`)
* `rolling.h`: rolling hash function (required by `BloomFilter.hpp` and `RollingHash.h`)
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
