libdir  := src/lib
bindir := src/bin
lib := $(wildcard $(libdir)/*.c)
CFLAGS := -Wall -Wextra -O3 

.PHONY: test
.PHONY: prof
.PHONY: prof_s
.PHONY: prof_t

build: $(lib)
	gcc $(CFLAGS) -fanalyzer $(bindir)/start.c -o start.o $(lib) -lgmp 
	gcc $(CFLAGS) -fanalyzer $(bindir)/test.c -o test.o $(lib) -lgmp -pthread
	clang $(CFLAGS) $(bindir)/start.c -o start.wasm $(lib) -lgmp 
	clang $(CFLAGS) $(bindir)/test.c -o test.wasm $(lib) -lgmp 
mm: $(lib)
	gcc $(CFLAGS) -fanalyzer -DMM_ENABLE $(bindir)/start.c -o start.o $(lib) -lgmp 
	gcc $(CFLAGS) -fanalyzer -DMM_ENABLE $(bindir)/test.c -o test.o $(lib) -lgmp 
	clang $(CFLAGS) -DMM_ENABLE $(bindir)/start.c -o start.wasm $(lib) -lgmp 
	clang $(CFLAGS) -DMM_ENABLE  $(bindir)/test.c -o test.wasm $(lib) -lgmp 
prof: $(lib)
	gcc -pg -Wall -Wextra -O2 $(bindir)/start.c -o prof.o $(lib) -lgmp 
	gcc -pg -Wall -Wextra -O2 $(bindir)/test.c -o prof_t.o $(lib) -lgmp 
test:
	./test.sh
bench:
	./bench.sh
clean:
	rm -f start.o test.o
	rm -f *.wasm
lint:
	bear -- make
format:
	./Cindent src/*.c
	./Cindent src/*.h
	rm src/*~
