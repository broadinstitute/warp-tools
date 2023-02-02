# The default target: the actual binaries the code is supposed to be built into.
all: bin/TagSort

# All tests produced by this Makefile. Add new tests you create to this list.
# UNIT TESTS MUST BE NAMED foo_test, with a corresponding foo_test.cpp in the
# test directory.
test: bin/alignment_datatype_test

CC = g++ -std=c++17 -Wall -Wno-sign-compare -O4

TAGSORT_OBJ = obj/alignment_datatype.o obj/metricgatherer.o obj/mitochondrial_gene_selector.o obj/partial_file_merge.o obj/partial_sort.o obj/tagsort_input_options.o

TAGSORT_LIBS = htslib/libhts.a -lz -lbz2 -llzma -lstdc++fs -lpthread

bin/TagSort: $(TAGSORT_OBJ) obj/tagsort_main.o
	$(CC) -o $@ $^ $(TAGSORT_LIBS)

obj/%.o: src/%.cpp src/*.h
	$(CC) -c -o $@ $<  -Ihtslib

.PHONY: clean
clean:
	rm -f obj/*.o bin/* *.o *.a

#=============GOOGLE TEST STUFF; DONT CHANGE===============
# Points to the root of Google Test, relative to where this file is.
# Remember to tweak this if you move this file.
GTEST_DIR = gtest

# Flags passed to the preprocessor.
# Set Google Test's header directory as a system directory, such that
# the compiler doesn't generate warnings in Google Test headers.
CPPFLAGS += -isystem $(GTEST_DIR)/include

# All Google Test headers.  Usually you shouldn't change this
# definition.
GTEST_HEADERS = $(GTEST_DIR)/include/gtest/*.h \
                $(GTEST_DIR)/include/gtest/internal/*.h

# Builds gtest.a and gtest_main.a.

# Usually you shouldn't tweak such internal variables, indicated by a
# trailing _.
GTEST_SRCS_ = $(GTEST_DIR)/src/*.cc $(GTEST_DIR)/src/*.h $(GTEST_HEADERS)

# For simplicity and to avoid depending on Google Test's
# implementation details, the dependencies specified below are
# conservative and not optimized.  This is fine as Google Test
# compiles fast and for ordinary users its source rarely changes.
gtest-all.o : $(GTEST_SRCS_)
	$(CC) $(CPPFLAGS) -I$(GTEST_DIR) -c \
            $(GTEST_DIR)/src/gtest-all.cc

gtest_main.o : $(GTEST_SRCS_)
	$(CC) $(CPPFLAGS) -I$(GTEST_DIR) -c \
            $(GTEST_DIR)/src/gtest_main.cc

gtest.a : gtest-all.o
	$(AR) $(ARFLAGS) $@ $^

gtest_main.a : gtest-all.o gtest_main.o
	$(AR) $(ARFLAGS) $@ $^

obj/%_test.o : test/%_test.cpp $(GTEST_HEADERS)
	$(CC) -Igtest/include -o $@ $< -c

bin/%_test : $(TAGSORT_OBJ) obj/%.o obj/%_test.o gtest_main.a
	$(CC) -lpthread $^ -o $@ $(TAGSORT_LIBS)

#==================END GOOGLE TEST STUFF===================
