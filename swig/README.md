# Swig

Make sure you have swig installed and included in your path.

To build a Perl5 module (run in swig/):
```
preinst-swig -Wall -c++ -perl5 BloomFilter.i
g++ -c BloomFilter_wrap.cxx -I/usr/lib64/perl5/CORE -fPIC -Dbool=char
g++ -Wall -shared BloomFilter_wrap.o -o BloomFilter.so
```

To run tests:
```
# To test insertions:
./writeBloom_rolling.pl
Usage: ./writeBloom_rolling.pl
-f  sequences to scaffold (Multi-FASTA format, required)
-k  k-mer value (default -k 15, optional)
-p  Bloom filter false positive rate (default -p 0.0001, optional - increase to prevent memory allocation errors)

# To test queries:
./testBloom_rolling.pl
Usage: ./testBloom_rolling.pl
-f  sequences to test (Multi-FASTA format, required)
-k  k-mer value (default -k 15, optional)
-n number of hashes (requred)
-s size of bloom filter in bits (reqired)
-b  Bloom filter
```

In order to compile, swig needs the following Perl5 headers:
```C++
#include "Extern.h"
#include "perl.h"
#include "XSUB.h"
```
If they are not located in /usr/lib64/perl5, you can run "perl -e 'use Config; print $Config{archlib};" to locate them.
