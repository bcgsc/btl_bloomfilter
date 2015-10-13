noinst_PROGRAMS = ReadProcessorTests \
	BloomFilterTests \
	BloomFilterInfoTests

ReadProcessorTests_LDADD = $(top_builddir)/BloomFilter/libbloomfilter.a -lz
ReadProcessorTests_SOURCES = ReadProcessorTests.cpp

BloomFilterTests_LDADD = $(top_builddir)/BloomFilter/libbloomfilter.a -lz
BloomFilterTests_SOURCES = BloomFilterTests.cpp

BloomFilterInfoTests_LDADD = $(top_builddir)/BloomFilter/libbloomfilter.a -lz
BloomFilterInfoTests_SOURCES = BloomFilterInfoTests.cpp