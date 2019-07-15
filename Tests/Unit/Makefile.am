# programs that will only be compiled with 'make check'
<<<<<<< HEAD
check_PROGRAMS = \
    BloomFilterTests \
    CountingBloomFilterTests
=======
check_PROGRAMS = BloomFilterTests BitVectorTests
>>>>>>> Makefile.am: add bit vector unit test to make file target

BloomFilterTests_SOURCES = BloomFilterTests.cpp
BloomFilterTests_CPPFLAGS = -I$(top_srcdir)
BloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

CountingBloomFilterTests_SOURCES = CountingBloomFilterTests.cpp
CountingBloomFilterTests_CPPFLAGS = -I$(top_srcdir)
CountingBloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

BitVectorTests_SOURCES = BloomFilterTests.cpp catch.hpp
BitVectorTests_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/vendor
BitVectorTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

# programs that will be run by 'make check' (after compiling)
TESTS = $(check_PROGRAMS)
