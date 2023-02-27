# The default target: the actual binaries the code is supposed to be built into.
all: bin/fastqprocess bin/fastq_slideseq bin/fastq_metrics bin/samplefastq

# All tests produced by this Makefile. Add new tests you create to this list.
# UNIT TESTS MUST BE NAMED foo_test, with a corresponding foo_test.cpp in the
# test directory.
test: bin/whitelist_corrector_test

CC = g++ -std=c++17 -Wall -Wno-sign-compare -O4

LIBS = -LlibStatGen -lStatGen -lz -lpthread -lstdc++fs -Lgzstream -lgzstream

COMMON_OBJ = obj/input_options.o obj/fastq_common.o obj/whitelist_corrector.o

bin/fastqprocess: $(COMMON_OBJ) obj/fastqprocess.o
	$(CC) -o $@ $^  $(LIBS)

bin/fastq_slideseq: $(COMMON_OBJ) obj/fastq_slideseq.o
	$(CC) -o $@ $^  $(LIBS)

bin/fastq_metrics: $(COMMON_OBJ) obj/fastq_metrics.o
	$(CC) -o $@ $^  $(LIBS)

bin/samplefastq: $(COMMON_OBJ) obj/samplefastq.o
	$(CC) -o $@ $^  $(LIBS)

obj/%.o: src/%.cpp src/*.h
	$(CC) -c -o $@ $<  -IlibStatGen/include -Igzstream

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

bin/%_test : $(COMMON_OBJ) obj/%.o obj/%_test.o gtest_main.a
	$(CC) -lpthread $^ -o $@ $(LIBS)

#==================END GOOGLE TEST STUFF===================
