CC       = gcc
CXX      = g++
LANGFLAG = -x c++
CPPFLAGS += -I slow5lib/include/
CFLAGS   += -g -Wall -O2 -std=c++11
LDFLAGS  += $(LIBS) -lpthread -lz -rdynamic
BUILD_DIR = build

ifeq ($(zstd),1)
LDFLAGS		+= -lzstd
endif

BINARY = sigfish
OBJ = $(BUILD_DIR)/main.o \
      $(BUILD_DIR)/dtw_main.o \
      $(BUILD_DIR)/sigfish.o \
      $(BUILD_DIR)/thread.o \
      $(BUILD_DIR)/events.o \
      $(BUILD_DIR)/model.o \
      $(BUILD_DIR)/cdtw.o \
	  $(BUILD_DIR)/genref.o \
	  $(BUILD_DIR)/cmain.o \
	  $(BUILD_DIR)/cfunc.o \
	  $(BUILD_DIR)/misc.o \
	  $(BUILD_DIR)/eval.o \

PREFIX = /usr/local
VERSION = `git describe --tags`

ifdef asan
	CFLAGS += -fsanitize=address -fno-omit-frame-pointer
	LDFLAGS += -fsanitize=address -fno-omit-frame-pointer
endif

.PHONY: clean distclean test

$(BINARY): $(OBJ) slow5lib/lib/libslow5.a
	$(CXX) $(CFLAGS) $(OBJ) slow5lib/lib/libslow5.a $(LDFLAGS) -o $@

$(BUILD_DIR)/main.o: src/main.c src/misc.h src/error.h src/sigfish.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/sigfish.o: src/sigfish.c src/misc.h src/error.h src/sigfish.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/thread.o: src/thread.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/dtw_main.o: src/dtw_main.c src/error.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/events.o: src/events.c src/misc.h src/ksort.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/model.o: src/model.c src/model.h  src/misc.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/cdtw.o: src/cdtw.c src/cdtw.h
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/genref.o: src/genref.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/cmain.o: src/cmain.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/cfunc.o: src/cfunc.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/misc.o: src/misc.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

$(BUILD_DIR)/eval.o: src/eval.c
	$(CXX) $(CFLAGS) $(CPPFLAGS) $(LANGFLAG) $< -c -o $@

slow5lib/lib/libslow5.a:
	$(MAKE) -C slow5lib zstd=$(zstd) no_simd=$(no_simd) zstd_local=$(zstd_local)  lib/libslow5.a

clean:
	rm -rf $(BINARY) $(BUILD_DIR)/*.o
	make -C slow5lib clean

# Delete all gitignored files (but not directories)
distclean: clean
	git clean -f -X
	rm -rf $(BUILD_DIR)/* autom4te.cache

test: $(BINARY)
	./sigfish dtw -g test/nCoV-2019.reference.fasta -s test/batch0.blow5  > test/res.paf
	diff -q test/batch0.exp.paf test/res.paf
	./sigfish dtw -g test/rnasequin_sequences_2.4.fa -s test/sequin_reads.blow5 --rna  > test/res_rna.paf
	diff -q test/res_rna.exp.paf test/res_rna.paf

eval: $(BINARY)
	./scripts/eval_dna.sh
	./scripts/eval_rna.sh


valgrind: $(BINARY)
	valgrind --leak-check=full ./sigfish dtw -g test/nCoV-2019.reference.fasta -s test/batch0.blow5 > test/res.paf
