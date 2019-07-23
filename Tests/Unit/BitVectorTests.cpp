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
#include "vendor/catch.hpp"

#include <assert.h>

using namespace std;
typedef uint64_t T;

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

	// Create a 2 bits vector
	const size_t filterSize = 1000000000;
	unsigned bitsPerCounter = 2;

	BitVector filter_2bits(filterSize, bitsPerCounter);

	/* Set up values 3,2,1 for element 0,1,LAST respectively */
	for (int i = 0; i < 3; i++) {
		filter_2bits.insert(0);
	}

	for (int i = 0; i < 2; i++) {
		filter_2bits.insert(1);
	}

	long lastElement = (filterSize * 8 / sizeof(T) / bitsPerCounter) - 1;
	filter_2bits.insert(lastElement);

	// Create a 4 bits vector
	bitsPerCounter = 4;

	BitVector filter_4bits(filterSize, bitsPerCounter);

	/* Set up values 15,2,1 for element 0,1,LAST respectively */
	for (int i = 0; i < 15; i++) {
		filter_4bits.insert(0);
	}

	for (int i = 0; i < 2; i++) {
		filter_4bits.insert(1);
	}

	lastElement = (filterSize * 8 / sizeof(T) / bitsPerCounter) - 1;
	filter_4bits.insert(lastElement);

	// Create a 8 bits vector
	bitsPerCounter = 8;

	BitVector filter_8bits(filterSize, bitsPerCounter);

	/* Set up values 255,2,1 for element 0,1,LAST respectively */
	for (int i = 0; i < 255; i++) {
		filter_8bits.insert(0);
	}

	for (int i = 0; i < 2; i++) {
		filter_8bits.insert(1);
	}

	lastElement = (filterSize * 8 / sizeof(T) / bitsPerCounter) - 1;
	filter_8bits.insert(lastElement);

	/* END COMMON SETUP CODE */

	SECTION("query elements")
	{
		/* Query elements 0,1,LAST and check for values 3,2,1, repectively */
		assert(filter_2bits[0] == 3);
		assert(filter_2bits[1] == 2);
		assert(filter_2bits[filter_2bits.size() - 1] == 1);

		/* Query elements 0,1,LAST and check for values 15,2,1, repectively */
		assert(filter_4bits[0] == 15);
		assert(filter_4bits[1] == 2);
		assert(filter_4bits[filter_4bits.size() - 1] == 1);

		/* Query elements 0,1,LAST and check for values 15,2,1, repectively */
		assert(filter_8bits[0] == 255);
		assert(filter_8bits[1] == 2);
		assert(filter_8bits[filter_8bits.size() - 1] == 1);
	}

	SECTION("query for no presence")
	{
		/* Query every elements besides 0,1,LAST and check they are equal to 0 */
		for (unsigned int i = 2; i < (filter_2bits.size() - 1); i++) {
			assert(filter_2bits[i] == 0);
		}

		/* Query every elements besides 0,1,LAST and check they are equal to 0 */
		for (unsigned int i = 2; i < (filter_4bits.size() - 1); i++) {
			assert(filter_4bits[i] == 0);
		}

		/* Query every elements besides 0,1,LAST and check they are equal to 0 */
		for (unsigned int i = 2; i < (filter_8bits.size() - 1); i++) {
			assert(filter_8bits[i] == 0);
		}
	}

	SECTION("insert elements and check")
	{
		/* Increment and check*/
		filter_2bits.insert(2);
		assert(filter_2bits[2] == 1);

		/* check that vector doesn't allow overflow */
		assert(filter_2bits.atomicIncrement(0) == false);
		/* check that vector value didn't cahnge */
		assert(filter_2bits[0] == 3);

		/* Increment and check*/
		filter_4bits.insert(2);
		assert(filter_4bits[2] == 1);

		/* check that vector doesn't allow overflow */
		assert(filter_4bits.atomicIncrement(0) == false);
		/* check that vector value didn't cahnge */
		assert(filter_4bits[0] == 15);

		/* Increment and check*/
		filter_8bits.insert(2);
		assert(filter_8bits[2] == 1);

		/* check that vector doesn't allow overflow */
		assert(filter_8bits.atomicIncrement(0) == false);
		/* check that vector value didn't cahnge */
		assert(filter_8bits[0] == 255);
	}
} /* end test fixture */
