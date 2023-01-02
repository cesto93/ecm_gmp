objdir  := src
bindir := src/bin
objects := $(addprefix $(objdir)/,m_ellcurv_fact.c m_ellcurv_lib.c mm_ellcurv_lib.c m_ellcurv_struct.c base_lib.c matbase_lib.c mpn_l.c)

.PHONY: test
.PHONY: prof
.PHONY: prof_s
.PHONY: prof_t

build: $(objects)
	 gcc -Wall -Wextra -O3 $(bindir)/start.c -o start.o $(objects) -lgmp -pthread
	 gcc -Wall -Wextra -O3 -static $(bindir)/test.c -o test.o $(objects) -lgmp -pthread
prof: $(objects)
	gcc -pg -Wall -Wextra -O2 $(bindir)/start.c -o prof.o $(objects) -lgmp -pthread
	gcc -pg -Wall -Wextra -O2 $(bindir)/test.c -o prof_t.o $(objects) -lgmp -pthread
test:
	./test.sh
bench:
	./bench.sh
clean:
	rm -f start.o test.o
lint:
	bear -- make
format:
	./Lindent src/*.c
	./Lindent src/*.h
