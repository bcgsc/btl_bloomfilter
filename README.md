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
