SUBDIRS = \
	Tests/AdHoc \
	Tests/Unit

include_HEADERS = \
	BloomFilter.hpp \
	CountingBloomFilter.hpp \
	KmerBloomFilter.hpp \
	MIBFQuerySupport.hpp \
	MIBloomFilter.hpp \
	nthash.hpp \
	ntHashIterator.hpp \
	ssHashIterator.hpp \
	stHashIterator.hpp

dist_doc_DATA = \
	README.md

EXTRA_DIST = autogen.sh

clang-format:
	for i in *.hpp *.h */*.cpp Tests/*/*.*pp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.hpp *.h */*.cpp Tests/*/*.*pp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi
