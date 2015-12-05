/*
 * BloomFilterTests.cpp
 * Unit Tests for hashmanager and bloomfilter classes
 *  Created on: Aug 14, 2012
 *      Author: cjustin
 */

/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN

/* lightweight unit test framework */
#include "catch.hpp"
#include "BloomFilter.hpp"

#include <string>
#include <assert.h>
#include <fstream>
#include <sstream>

using namespace std;

TEST_CASE("test fixture", "[BloomFilter]")
{
	/*
	 * NOTES:
	 * - The SECTION blocks below are separate tests that share the
	 * same setup code
	 * - The common setup code is _re-run_ before each SECTION
	 * - In unit test terminology, this type of setup is known as a
	 * "test fixture"
	 * - See https://github.com/philsquared/Catch/blob/master/docs/tutorial.md#test-cases-and-sections
	 * for details
	 */

	/* START COMMON SETUP CODE */

	//test Bloom filter with some elements
	size_t filterSize = 1000000000;
	BloomFilter filter(filterSize, 5, 20);
	filter.insert("ATCGGGTCATCAACCAATAT");
	filter.insert("ATCGGGTCATCAACCAATAC");
	filter.insert("ATCGGGTCATCAACCAATAG");
	filter.insert("ATCGGGTCATCAACCAATAA");

	/* END COMMON SETUP CODE */

	SECTION("query elements")
	{
		//Check if filter is able to report expected results
		REQUIRE(filter.contains("ATCGGGTCATCAACCAATAT"));
		REQUIRE(filter.contains("ATCGGGTCATCAACCAATAC"));
		REQUIRE(filter.contains("ATCGGGTCATCAACCAATAG"));
		REQUIRE(filter.contains("ATCGGGTCATCAACCAATAA"));

		REQUIRE(!filter.contains("ATCGGGTCATCAACCAATTA"));
		REQUIRE(!filter.contains("ATCGGGTCATCAACCAATTC"));
	}

	SECTION("save/load Bloom file")
	{
		//Check storage can occur properly
		string filename = "/tmp/bloomFilter.bf";
		filter.storeFilter(filename);
		ifstream ifile(filename.c_str());
		assert(ifile.is_open());
		ifile.seekg(0, ios::end); // move to end of file
		size_t fileSize = ifile.tellg(); // file size in bytes
		//file size should be same as filter size (Round to block size)
		if (filterSize % 64 > 0) {
			assert((filterSize + (64 - (filterSize% 64))) == fileSize*8);
		} else {
			assert(filterSize == fileSize*8);
		}
		ifile.close();

		//check loading of stored filter
		BloomFilter filter2(filterSize, 5, 20, filename);

		//Check if loaded filter is able to report expected results
		REQUIRE(filter2.contains("ATCGGGTCATCAACCAATAT"));
		REQUIRE(filter2.contains("ATCGGGTCATCAACCAATAC"));
		REQUIRE(filter2.contains("ATCGGGTCATCAACCAATAG"));
		REQUIRE(filter2.contains("ATCGGGTCATCAACCAATAA"));

		REQUIRE(!filter2.contains("ATCGGGTCATCAACCAATTA"));
		REQUIRE(!filter2.contains("ATCGGGTCATCAACCAATTC"));

		remove(filename.c_str());
	}

} /* end test fixture */

