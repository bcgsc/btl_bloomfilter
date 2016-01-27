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
./test.pl
```

In order to compile, swig needs the following Perl5 headers:
```C++
#include "Extern.h"
#include "perl.h"
#include "XSUB.h"
```
If they are not located in /usr/lib64/perl5, you can run "perl -e 'use Config; print $Config{archlib};" to locate them.
