# programs that will not be installed with 'make install'
noinst_PROGRAMS = BloomFilterTests \
	BloomFilterAdHocTests \
	BloomFilterInfoTests \
	ParallelFilter

# programs that will only be compiled with 'make check'
check_PROGRAMS = BloomFilterTests \
	BloomFilterAdHocTests \
	BloomFilterInfoTests \
	ParallelFilter

BloomFilterTests_SOURCES = BloomFilterTests.cpp
BloomFilterTests_CPPFLAGS = -I$(top_srcdir)

BloomFilterAdHocTests_SOURCES = BloomFilterAdHocTests.cpp
BloomFilterAdHocTests_CPPFLAGS = -I$(top_srcdir)
BloomFilterAdHocTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

BloomFilterInfoTests_SOURCES = BloomFilterInfoTests.cpp
BloomFilterInfoTests_CPPFLAGS = -I$(top_srcdir)

ParallelFilter_SOURCES = ParallelFilter.cpp
ParallelFilter_CPPFLAGS = -I$(top_srcdir)
ParallelFilter_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

# programs that will be run by 'make check' (after compiling)
TESTS = BloomFilterTests
