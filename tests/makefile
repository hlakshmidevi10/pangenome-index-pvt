# Path to sdsl-lite
SDSL_DIR ?= /Users/seeskand/Documents/sdsl-lite
include $(SDSL_DIR)/Make.helper

BUILD_BIN=bin
BUILD_LIB=lib
BUILD_OBJ=obj
SOURCE_DIR=src

# Compiler - default to g++ for Linux, clang++ for macOS
MY_CXX ?= g++

# Parallelization flags
PARALLEL_FLAGS=-fopenmp -pthread

# Libraries
LIBS=-L$(LIB_DIR) -L../deps/grlBWT/build -lgbwt -lgbwtgraph -lhandlegraph -lsdsl -ldivsufsort -ldivsufsort64 -lgrlbwt -lgtest -lgtest_main

# macOS-specific OpenMP & compiler handling
ifeq ($(shell uname -s), Darwin)
    MY_CXX := clang++
    HOMEBREW_PREFIX := /opt/homebrew

    # The compiler complains about -fopenmp instead of missing input.
    ifeq ($(strip $(shell $(MY_CXX) -fopenmp /dev/null -o/dev/null 2>&1 | grep fopenmp | wc -l)), 1)
        $(info Detected Apple Clang; configuring OpenMP with libomp...)

        # The compiler only needs to do the preprocessing.
        PARALLEL_FLAGS = -Xpreprocessor -fopenmp -pthread

        # Find libomp installed by Homebrew or MacPorts.
        ifeq ($(shell if [ -e $(HOMEBREW_PREFIX)/include/omp.h ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp via Homebrew.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/include
            LIBS += -L$(HOMEBREW_PREFIX)/lib
        else ifeq ($(shell if [ -d $(HOMEBREW_PREFIX)/opt/libomp/include ]; then echo 1; else echo 0; fi), 1)
            $(info Found keg-only libomp via Homebrew.)
            PARALLEL_FLAGS += -I$(HOMEBREW_PREFIX)/opt/libomp/include
            LIBS += -L$(HOMEBREW_PREFIX)/opt/libomp/lib
        else ifeq ($(shell if [ -d /opt/local/lib/libomp ]; then echo 1; else echo 0; fi), 1)
            $(info Found libomp via MacPorts.)
            PARALLEL_FLAGS += -I/opt/local/include/libomp
            LIBS += -L/opt/local/lib/libomp
        else
            $(error Could not find libomp. Please install it using Homebrew or MacPorts.)
        endif

        # We also need to link it.
        LIBS += -lomp
    endif

    CXX_FLAGS += -I$(HOMEBREW_PREFIX)/include
endif

CXX_FLAGS=$(MY_CXX_FLAGS) $(PARALLEL_FLAGS) $(MY_CXX_OPT_FLAGS) -I../include -I$(INC_DIR) -I../lib -I../deps/r-index/internal -I../deps/vg -I../deps/grlBWT/include

HEADERS=$(wildcard ../include/pangenome_index/*.hpp)

#LIBOBJS=$(addprefix $(BUILD_OBJ)/,tag_arrays.o r-index.o)
#LIBRARY=$(BUILD_LIB)/libpanindexer.a

PROGRAMS=test_rindex test_bplus_tree benchmark_query_tags

LIBRARY=../lib/libpanindexer.a

.PHONY: all clean test
all:$(PROGRAMS)

directories: $(BUILD_BIN) $(BUILD_LIB) $(BUILD_OBJ)


%.o:%.cpp $(HEADERS)
	$(MY_CXX) $(CXX_FLAGS) -c $<

$(PROGRAMS):%:%.o $(LIBRARY)
	$(MY_CXX) $(CXX_FLAGS) -o $@ $< $(LIBRARY) $(LIBS)


test: $(PROGRAMS)
	$(foreach PROG,$(PROGRAMS),./$(PROG) || exit 1;)

clean:
	rm -f $(PROGRAMS) *.o
