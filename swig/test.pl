#!/usr/bin/perl

use BloomFilter;

#Same inputs as the adhoc tests
$filterSize = 1000000000;

#Check if filter is able to report expected results
$filter = new BloomFilter::BloomFilter($filterSize, 5, 20);

BloomFilter::BloomFilter::insert($filter, "ATCGGGTCATCAACCAATAT");
BloomFilter::BloomFilter::insert($filter, "ATCGGGTCATCAACCAATAC");
BloomFilter::BloomFilter::insert($filter, "ATCGGGTCATCAACCAATAG");
BloomFilter::BloomFilter::insert($filter, "ATCGGGTCATCAACCAATAA");

if (!BloomFilter::BloomFilter::contains($filter, "ATCGGGTCATCAACCAATAT")
	&&!BloomFilter::BloomFilter::contains($filter, "ATCGGGTCATCAACCAATAC")
	&& !BloomFilter::BloomFilter::contains($filter, "ATCGGGTCATCAACCAATAG")
	&& !BloomFilter::BloomFilter::contains($filter, "ATCGGGTCATCAACCAATAA")) {
	print "Filter did not contain expected. \n";
}

if (BloomFilter::BloomFilter::contains($filter, "ATCGGGTCATCAACCAATTA")
	&& BloomFilter::BloomFilter::contains($filter, "ATCGGGTCATCAACCAATTC")) {
	print "Filter contained unexpected. \n";
}

print "de novo bf tests done \n";

#Check storage can occur properly
$fileName = "BloomFilter.bf";
BloomFilter::BloomFilter::storeFilter($filter, $fileName);

$filter2 = new BloomFilter::BloomFilter($filterSize, 5, 20, $fileName);

if (!BloomFilter::BloomFilter::contains($filter2, "ATCGGGTCATCAACCAATAT")
	&&!BloomFilter::BloomFilter::contains($filter2, "ATCGGGTCATCAACCAATAC")
	&& !BloomFilter::BloomFilter::contains($filter2, "ATCGGGTCATCAACCAATAG")
	&& !BloomFilter::BloomFilter::contains($filter2, "ATCGGGTCATCAACCAATAA")) {
	print "Filter2 did not contain expected. \n";
}

if (BloomFilter::BloomFilter::contains($filter2, "ATCGGGTCATCAACCAATTA")
	&& BloomFilter::BloomFilter::contains($filter2, "ATCGGGTCATCAACCAATTC")) {
	print "Filter2 contained unexpected. \n";

}
print "premade bf tests done\n";

#Tests to see if precomputed inserts and contains can be called
$a = new BloomFilter::SizetVector(5);
BloomFilter::SizetVector::push($a, 1);
BloomFilter::SizetVector::push($a, 2);
BloomFilter::SizetVector::push($a, 3);
BloomFilter::SizetVector::push($a, 4);
BloomFilter::SizetVector::push($a, 5);

BloomFilter::BloomFilter::insert($filter, $a);
BloomFilter::BloomFilter::contains($filter, $a);

print "Done!\n";
exit;
