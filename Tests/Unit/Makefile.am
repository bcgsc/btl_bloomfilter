# programs that will only be compiled with 'make check'
check_PROGRAMS = \
    BitVectorTests \
    BloomFilterTests \
    CountingBloomFilterTests

BitVectorTests_SOURCES = BloomFilterTests.cpp catch.hpp
BitVectorTests_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/vendor
BitVectorTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

BloomFilterTests_SOURCES = BloomFilterTests.cpp catch.hpp
BloomFilterTests_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/vendor
BloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

CountingBloomFilterTests_SOURCES = CountingBloomFilterTests.cpp catch.hpp
CountingBloomFilterTests_CPPFLAGS = -I$(top_srcdir) -I$(top_srcdir)/vendor
CountingBloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS) -std=c++11

# programs that will be run by 'make check' (after compiling)
TESTS = $(check_PROGRAMS)
