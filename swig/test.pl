#!/usr/bin/perl

use BloomFilter;

#Same inputs as the adhoc tests
$filterSize = 1000000000;

#Check if filter is able to report expected results
$filter = BloomFilter::BloomFilter->new($filterSize, 5, 20);

$filter->insert("ATCGGGTCATCAACCAATAT");
$filter->insert("ATCGGGTCATCAACCAATAC");
$filter->insert("ATCGGGTCATCAACCAATAG");
$filter->insert("ATCGGGTCATCAACCAATAA");

if (!$filter->contains("ATCGGGTCATCAACCAATAT")
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

#RollingHashIterator tests
my $k = 5;
$str = "TAGAATCACCCAAAGA";
$bloom = new BloomFilter::BloomFilter(10000, 4, $k);
$itr = new BloomFilter::RollingHashIterator($str, 4, $k);

BloomFilter::insertSeq($bloom, $str, 4, $k);

#my $count = 0;
#my $next = $itr->getNext();
#while ($next) {
#	$bloom->insert($next);
#    print substr($str, $count, $k) . " " .  $count++ . "\n";
#	$next =  $itr->getNext();
#}

for (my $i = 0; $i < length($str) - $k + 1; $i++) {
	my $kmer = substr($str, $i, $k);
	print $i . " ";
	if($bloom->contains($kmer)){
		print $kmer . " found\n";
	}
	else{
		print $kmer . " not found\n";
	}
}

print "Done!\n";

exit;
