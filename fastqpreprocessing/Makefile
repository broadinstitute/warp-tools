CC = g++ -std=c++17 -Wall -Wno-sign-compare -O4

LIBS = -LlibStatGen -lStatGen -lz -lpthread -lstdc++fs -Lgzstream -lgzstream

install: bin/fastqprocess bin/TagSort bin/fastq_slideseq bin/fastq_metrics bin/samplefastq

COMMON_OBJ = obj/utilities.o obj/input_options.o obj/fastq_common.o

TAGSORT_OBJ = obj/tagsort.o obj/htslib_tagsort.o obj/metricgatherer.o obj/input_options.o

obj/%.o: src/%.cpp src/*.h
	$(CC) -c -o $@ $<  -IlibStatGen/include -Ihtslib -Igzstream

bin/fastqprocess: $(COMMON_OBJ) obj/fastqprocess.o
	$(CC) -o $@ $^  $(LIBS)

bin/TagSort: $(TAGSORT_OBJ)
	$(CC) -o $@ $(TAGSORT_OBJ) htslib/libhts.a -lz -lbz2 -llzma -Lgzstream -lgzstream -lstdc++fs -lpthread

bin/fastq_slideseq: $(COMMON_OBJ) obj/fastq_slideseq.o
	$(CC) -o $@ $^  $(LIBS)

bin/fastq_metrics: $(COMMON_OBJ) obj/fastq_metrics.o
	$(CC) -o $@ $^  $(LIBS)

bin/samplefastq: $(COMMON_OBJ) obj/samplefastq.o
	$(CC) -o $@ $^  $(LIBS)

.PHONY: clean
clean:
	rm -f obj/*.o bin/*
