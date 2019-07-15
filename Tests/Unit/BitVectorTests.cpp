/*
 * BitVector.cpp
 * Unit Tests for BitVector class
 *  Created on: July 15, 2019
 *      Author: Johnathan Wong
 */

/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN

/* lightweight unit test framework */
#include "BitVector.hpp"
#include "catch.hpp"

#include <assert.h>

using namespace std;

TEST_CASE("test fixture", "[BitVector]")
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
	const unsigned bitsPerCounter = 2;

	BitVector<uint8_t> filter(filterSize, bitsPerCounter);

	/* Set up values 3,2,1 for element 0,1,LAST respectively */

	for (int i = 0; i < 3; i++) {
		filter.atomicIncrement(0);
	}

	for (int i = 0; i < 2; i++) {
		filter.atomicIncrement(1);
	}

	long lastElement = (filterSize * 8 * sizeof(uint8_t)/ 2 ) - 1 ;
	for (int i = 0; i < 1; i++) {
		filter.atomicIncrement(lastElement);
	}

	/* END COMMON SETUP CODE */

	SECTION("query elements")
	{
		/* Query elements 0,1,LAST for values 3,2,1, repectively */
        assert(filter[0] == 3);
        assert(filter[1] == 2);
        assert(filter[filter.size() - 1] == 1);
	}

	SECTION("query for no presence")
	{
		/* Query everyelements besides 0,1,LAST and check they are equal to 0 */
        for (int i = 2; i < (filter.size() - 1); i++) {
            assert(filter[i] == 0);
        }
	}

	SECTION("insert elements and check")
	{
		/* Increment and check*/
		filter.atomicIncrement(2);
        assert(filter[2] == 1);

        /* check that vector doesn't allow overflow */
        assert(filter.atomicIncrement(0) == false);
        /* check that vector value didn't cahnge */
        assert(filter[0] == 3);
	}
} /* end test fixture */
