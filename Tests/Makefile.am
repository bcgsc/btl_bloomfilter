check_PROGRAMS = ReadProcessorTests \
	BloomFilterTests \
	BloomFilterCategorizerTests \
	BloomFilterMakerTests \
	WindowedParser \
	BloomFilterInfoTests

ReadProcessorTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Bloomfilter/libbloomfilter.a -lz
ReadProcessorTests_SOURCES = ReadProcessorTests.cpp

BloomFilterTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Bloomfilter/libbloomfilter.a -lz
BloomFilterTests_SOURCES = BloomFilterTests.cpp

BloomFilterInfoTests_LDADD = $(top_builddir)/DataLayer/libdatalayer.a \
	$(top_builddir)/Bloomfilter/libbloomfilter.a -lz
BloomFilterInfoTests_SOURCES = BloomFilterInfoTests.cpp