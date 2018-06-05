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
#include "ntHashIterator.hpp"

#include <string>
#include <assert.h>
#include <fstream>
#include <sstream>
#include <cstring>
#include <stdlib.h>

using namespace std;

/** create a uniquely-named temp file (tedious!) */
string createTempFile()
{
	const unsigned MAX_FILENAME_SIZE = 1024;
	const char* fileTemplate = "/XXXXXX";
	char filename[MAX_FILENAME_SIZE + 1];

	/* allow override of default tmp dir */
	char* tmpdir = getenv("TMPDIR");
	if (tmpdir)
		strcpy(filename, tmpdir);
	else
		strcpy(filename, "/tmp");

	assert(strlen(filename) + strlen(fileTemplate) <= MAX_FILENAME_SIZE);
	strcat(filename, fileTemplate);

	int fd = mkstemp(filename);
	if (fd == -1) {
		perror("failed to create temp file");
		exit(EXIT_FAILURE);
	}
	close(fd);

	return string(filename);
}

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

	const size_t filterSize = 1000000000;
	const unsigned numHashes = 5;
	const unsigned k = 4;
	const char* seq = "ACGTAC";

	BloomFilter filter(filterSize, numHashes, k);

	/* insert k-mers ACGT, CGTA, GTAC */

	ntHashIterator insertIt(seq, numHashes, k);
	while(insertIt != insertIt.end()) {
		filter.insert(*insertIt);
		++insertIt;
	}

	/* END COMMON SETUP CODE */

	SECTION("query elements")
	{
		/* check that k-mers were correctly inserted */

		ntHashIterator queryIt(seq, numHashes, k);
		while(queryIt != queryIt.end()) {
			assert(filter.contains(*queryIt));
			++queryIt;
		}
	}

	SECTION("save/load Bloom file")
	{
		/* write filter */

		string filename = createTempFile();
		filter.storeFilter(filename);
		ifstream ifile(filename.c_str());

		/* check size of newly-created file */

		assert(ifile.is_open());
		ifile.seekg(0, ios::end); // move to end of file
		size_t fileSize = ifile.tellg(); // file size in bytes
		//file size should be same as filter size (Round to block size)
		if (filterSize % 64 > 0) {
			assert((filterSize + (64 - (filterSize% 64))) + sizeof(BloomFilter::FileHeader)*8 == fileSize*8);
		} else {
			assert(filterSize + sizeof(BloomFilter::FileHeader)*8 == fileSize*8);
		}
		ifile.close();

		/* check loading of stored filter */

		BloomFilter filter2(filename);

		/* check if loaded filter is able to report expected results */

		ntHashIterator queryIt(seq, numHashes, k);
		while(queryIt != queryIt.end()) {
			assert(filter.contains(*queryIt));
			++queryIt;
		}

		/* cleanup */

		remove(filename.c_str());
	}

} /* end test fixture */

