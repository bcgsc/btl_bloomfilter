noinst_PROGRAMS = BloomFilterTests \
	BloomFilterInfoTests \
	ParallelFilter

BloomFilterTests_SOURCES = BloomFilterTests.cpp
BloomFilterTests_CPPFLAGS = -I$(top_srcdir)
BloomFilterTests_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)

BloomFilterInfoTests_SOURCES = BloomFilterInfoTests.cpp
BloomFilterInfoTests_CPPFLAGS = -I$(top_srcdir)

ParallelFilter_SOURCES = ParallelFilter.cpp
ParallelFilter_CPPFLAGS = -I$(top_srcdir)
ParallelFilter_CXXFLAGS = $(AM_CXXFLAGS) $(OPENMP_CXXFLAGS)
