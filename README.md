# bloomfilter
The BTL C/C++ Common Bloom filters for bioinformatics projects, as well as any APIs created for other programming languages.

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

* magic
	Description: bf magic string 
	Type: char[8]
	Value: BLOOMFXX
* hlen
	Description: length of the header text
	Type: uint32_t
	Value:
* header
	Description: Plain header text
	Type: char[hlen]
	Value:
	- size:
		Description: The size of Bloom filter
		Type: uint64_t
		Value:
	- hmethod
		Description: hashing algorithm
		Type: char[8]
		Value:
	- nhash
		Description: number of hashes
		Type: uint32_t
		Value:
	- kmer [optional]
		Description: k-mer size
		Type: uint32_t
		Value:
	- seed [optional]
		Description: initial seeds for different hashes
		Type: uint64_t[nhash]
		Value: [0,1, ..., nhash-1] 
	- dFPR [optional]
		Description: desired false positve rate
		Type: double
		Value:
	- aFPR [optional]
		Description: approximate false positive rate
		Type: double
		Value:
	- rFPR [optional]
		Description: redundant false positive rate
		Type: double
		Value:
	- nEntry [optional]
		Description: number of entries
		Type: uint64_t
		Value:
	- tEntry [optional]
		Description: total number of entries
		Type: uint64_t
		Value:
* filter 
	Description: Bloom filter content
	Type: uchar[size]
