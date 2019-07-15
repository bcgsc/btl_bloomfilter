/*
 * CountingBloomFilterTests.cpp
 * Unit Tests for CountingBloomFilter class
 *  Created on: July 15, 2019
 *      Author: Johnathan Wong
 */

/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN

/* lightweight unit test framework */
#include "CountingBloomFilter.hpp"
#include "catch.hpp"
#include "ntHashIterator.hpp"

#include <assert.h>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <string>

using namespace std;

/** create a uniquely-named temp file (tedious!) */
string
createTempFile()
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

TEST_CASE("test fixture", "[CountingBloomFilter]")
{
	/*
	 * NOTES:
	 * - The SECTION blocks below are separate tests that share the
	 * same setup code
	 * - The common setup code is _re-run_ before each SECTION
	 * - In unit test terminology, this type of setup is known as a
	 * "test fixture"
	 * - See
	 * https://github.com/philsquared/Catch/blob/master/docs/tutorial.md#test-cases-and-sections for
	 * details
	 */

	/* START COMMON SETUP CODE */

	const size_t filterSize = 1000000000;
	const unsigned numHashes = 5;
	const unsigned threshold = 1;
	const unsigned k = 4;
	const char* seq = "ACGTAC";

	CountingBloomFilter<uint8_t> filter(filterSize, numHashes, k, threshold);

	/* insert k-mers ACGT, CGTA, GTAC */

	ntHashIterator insertIt(seq, numHashes, k);
	while (insertIt != insertIt.end()) {
		filter.insert(*insertIt);
		++insertIt;
	}

	/* END COMMON SETUP CODE */

	SECTION("query elements")
	{
		/* check that k-mers were correctly inserted */

		ntHashIterator queryIt(seq, numHashes, k);
		while (queryIt != queryIt.end()) {
			assert(filter.contains(*queryIt));
			++queryIt;
		}
	}

	SECTION("save/load Bloom file")
	{
		/* write filter */

		string filename = createTempFile();
		filter.writeFilter(filename);
		ifstream ifile(filename.c_str());

		/* check size of newly-created file */

		assert(ifile.is_open());
        /* File header has no fixed element but "[HeaderEnd]"
           so check for "[HeaderEnd]" in file*/
        std::string headerEnd = "[HeaderEnd]";
        std::string line;
        bool headerEndCheck = false;
        while (std::getline(ifile, line)) {
            if (line == headerEnd) {
                headerEndCheck = true;
                break;
            }
        }
        assert(headerEndCheck);

        size_t currPos = ifile.tellg();
        ifile.seekg(0, ios::end);
        size_t endPos = ifile.tellg();
		// file size - header should be same as filter size
        assert((endPos - currPos) == filterSize);
		ifile.close();

		/* check loading of stored filter */

		CountingBloomFilter<uint8_t> filter2(filename, threshold);

		/* check if loaded filter is able to report expected results */

		ntHashIterator queryIt(seq, numHashes, k);
		while (queryIt != queryIt.end()) {
			assert(filter2.contains(*queryIt));
			++queryIt;
		}

		/* cleanup */

		remove(filename.c_str());
	}

} /* end test fixture */
