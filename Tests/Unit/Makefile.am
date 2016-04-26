# programs that will only be compiled with 'make check'
check_PROGRAMS = BloomFilterTests BloomMapTests CountingBloomFilterTests

BloomFilterTests_SOURCES = BloomFilterTests.cpp
BloomFilterTests_CPPFLAGS = -I$(top_srcdir)
BloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

BloomMapTests_SOURCES = BloomMapTests.cpp
BloomMapTests_CPPFLAGS = -I$(top_srcdir)
BloomMapTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

CountingBloomFilterTests_SOURCES = CountingBloomFilterTests.cpp
CountingBloomFilterTests_CPPFLAGS = -I$(top_srcdir)
CountingBloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

# programs that will be run by 'make check' (after compiling)
TESTS = $(check_PROGRAMS)
