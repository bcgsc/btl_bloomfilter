/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN
/* lightweight unit test framework */
#include "catch.hpp"
#include <vector>

using namespace std;

/* remove #if/#endif when CountingBloomFilter is implemented */
#include "CountingBloomFilter.hpp"

TEST_CASE("insert and query", "[CountingBloomFilter]")
{
	const size_t bloomSize = 1000;
	const unsigned numHashes = 4;

	CountingBloomFilter<int> countingBloom(bloomSize, numHashes);

	vector<size_t> hashes;
	vector<size_t> elements;

	hashes.push_back(0);
	hashes.push_back(1);
	hashes.push_back(2);
	hashes.push_back(3);

	REQUIRE(countingBloom[hashes] == 0);
	countingBloom.insert(hashes);
	REQUIRE(countingBloom[hashes] == 1);
	countingBloom.insert(hashes);
	REQUIRE(countingBloom[hashes] == 2);
	countingBloom.insert(hashes);
	REQUIRE(countingBloom[hashes] == 3);

}
