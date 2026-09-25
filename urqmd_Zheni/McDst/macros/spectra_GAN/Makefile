# Define C++ compiler
CXX = g++

# Define linker
LINKER = g++

# Debug flag: either nothing or -ggdb
DEBUGFLAG = -ggdb

# Opitmization flag
OPTFLAG = -O2

# Path to McDst (either absolute or relative)
MCDSTPATH = /star/u/annakraeva/urqmd_hbt_Zheni/urqmd_Zheni/McDst

# ROOT flag: -D__ROOT__
ROOTFLAG =

# Define C++ flags
CXXFLAGS = $(shell root-config --cflags) -I. -I$(MCDSTPATH) $(OPTFLAG) -fPIC -Wall -pipe -std=c++17
CXXFLAGS += $(DEBUGFLAG)

# Define libraries to be added
LIBS = $(shell root-config --libs) -lEG
LIBS += -L$(MCDSTPATH) -lMcDst

# Define path to include directories
INCS = $(shell root-config --incdir)

# Linker flags
LDFLAGS = $(shell root-config --cflags) -I. -I$(MCDSTPATH) $(OPTFLAG) -fPIC -Wall -pipe

all: spectraFromMcDst 

# Compile C++ files
spectraFromMcDst: spectraFromMcDst.o
	$(LINKER) $< $(LDFLAGS) $(LIBS) $(GLIBS) -o $@

# Compile C++ files
%.o: %.C
	$(CXX) $(CXXFLAGS) -c -o $@ $<
	@echo 'Compiling: $@ $<'

# PHONY options
.PHONY: clean distclean

# Remove object files and dictionary
clean:
	rm -vf *.o

# Remove object files, dictionary and library
distclean:
	rm -vf *.o spectraFromMcDst 

