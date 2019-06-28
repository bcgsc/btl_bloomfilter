SUBDIRS = \
	Tests/AdHoc \
	Tests/Unit

include_HEADERS = \
	BloomFilter.hpp \
	CountingBloomFilter.hpp \
	KmerBloomFilter.hpp \
	MIBFQuerySupport.hpp \
	MIBloomFilter.hpp \
	vendor/nthash.hpp \
	vendor/ntHashIterator.hpp \
	vendor/ssHashIterator.hpp \
	vendor/stHashIterator.hpp

dist_doc_DATA = \
	README.md

EXTRA_DIST = autogen.sh


# Check the C++ source code for white-space errors with clang-format.

clang-format:
	for i in *.hpp *.h */*.cpp Tests/*/*.*pp; do clang-format -style=file $$i >$$i.fixed; done
	for i in *.hpp *.h */*.cpp Tests/*/*.*pp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi

# Check the C++ source code for errors with clang-tidy.

clang-tidy:
	clang-tidy -warnings-as-errors='*' *.hpp *.h */*.cpp Tests/*/*.*pp -- -std=c++11 -x c++ -I./vendor