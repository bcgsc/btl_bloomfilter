# Swig

To compile perl module (run in swig/ dir):

/home/rwarren/bin/swig-3.0.7/preinst-swig -Wall -c++ -perl5 BloomFilter.i

g++ -c BloomFilter_wrap.cxx -I/usr/lib64/perl5/CORE -fPIC -Dbool=char

g++ -Wall -shared BloomFilter_wrap.o -o BloomFilter.so

To run tests:

./test.pl
