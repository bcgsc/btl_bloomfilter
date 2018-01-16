/* automatically create main() function to run tests */
#define CATCH_CONFIG_MAIN
/* lightweight unit test framework */
#include "catch.hpp"
#include <sdsl/int_vector.hpp>

using namespace std;

#include "MIBloomFilter.hpp"

TEST_CASE("TEST EXAMPLE", "[MIBloomFilter]")
{
	sdsl::bit_vector bv(1000);
	cout << bv.size() << end;
}
