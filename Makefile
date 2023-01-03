libdir  := src/lib
bindir := src/bin
lib := $(addprefix $(libdir)/,m_ellcurv.c mm_ellcurv.c m_ellcurv_struct.c base_lib.c matbase_lib.c mpn_l.c)

.PHONY: test
.PHONY: prof
.PHONY: prof_s
.PHONY: prof_t

build: $(lib)
	gcc -Wall -Wextra -O3 $(bindir)/start.c -o start.o $(lib) -lgmp -pthread
	gcc -Wall -Wextra -O3 -static $(bindir)/test.c -o test.o $(lib) -lgmp -pthread
	clang $(bindir)/start.c -o start.wasm $(lib) -lgmp -pthread
	clang $(bindir)/test.c -o test.wasm $(lib) -lgmp -pthread
prof: $(lib)
	gcc -pg -Wall -Wextra -O2 $(bindir)/start.c -o prof.o $(lib) -lgmp -pthread
	gcc -pg -Wall -Wextra -O2 $(bindir)/test.c -o prof_t.o $(lib) -lgmp -pthread
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
