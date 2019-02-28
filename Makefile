#Set inclusion paths here (if boost, bamtools, or args are installed outside your path)
INCLUDE_DIRS=-Irnaseqc -Irnaseqc/src -Irnaseqc/SeqLib -Irnaseqc/SeqLib/htslib/
#Set library paths here (if boost or bamtools are installed outside your path)
LIBRARY_PATHS=
#Set to 0 if you encounter linker errors regarding strings from the bamtools library
ABI=1
#Provide full paths here to .a archives for libraries which should be statically linked
STATIC_LIBS=
#List of remaining libraries that will be dynamically linked
LIBS=-lboost_filesystem -lboost_regex -lboost_system -lz -llzma -lbz2 -lpthread

CC=g++
STDLIB=-std=c++14
CFLAGS=-Wall $(STDLIB) -D_GLIBCXX_USE_CXX11_ABI=$(ABI) -O3
SOURCES=scrinvex.cpp
SRCDIR=src
OBJECTS=$(SOURCES:.cpp=.o)
SEQFLAGS=$(STDLIB) -D_GLIBCXX_USE_CXX11_ABI=$(ABI)

scrinvex: $(foreach file,$(OBJECTS),$(SRCDIR)/$(file)) rnaseqc/rnaseqc.a rnaseqc/SeqLib/bin/libseqlib.a rnaseqc/SeqLib/bin/libhts.a
	$(CC) -O3 $(LIBRARY_PATHS) -o $@ $^ $(STATIC_LIBS) $(LIBS)

%.o: %.cpp
	$(CC) $(CFLAGS) -I. $(INCLUDE_DIRS) -c -o $@ $<

rnaseqc/SeqLib/bin/libseqlib.a rnaseqc/SeqLib/bin/libhts.a:
	cd rnaseqc/SeqLib && ./configure && make CXXFLAGS="$(SEQFLAGS)" && make install

rnaseqc/rnaseqc.a:
	cd rnaseqc && make lib ABI=$(ABI)

.PHONY: clean

clean:
	rm $(wildcard $(SRCDIR)/*.o) || echo "Nothing to clean in scrinvex"
	cd rnaseqc && make clean || echo "Nothing to clean in RNA-SeQC"
	cd rnaseqc/SeqLib && make clean || echo "Nothing to clean in SeqLib"
