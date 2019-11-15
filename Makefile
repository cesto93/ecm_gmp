objdir  := source
objects := $(addprefix $(objdir)/,m_ellcurv_fact.c m_ellcurv_lib.c mm_ellcurv_lib.c m_ellcurv_struct.c base_lib.c matbase_lib.c mpn_l.c)

.PHONY: test
.PHONY: prof
.PHONY: prof_s
.PHONY: prof_t

install_mac: $(objects)
	 gcc -Wall -Wextra -O3 source/start.c -o start.o $(objects) -lgmp -pthread
install: $(objects)
	 gcc -Wall -Wextra -O3 -static source/start.c -o start.o $(objects) -lgmp -pthread
debug_s: $(objects)
	gcc -g -Wall -Wextra -O2 source/start.c -o debug.o $(objects) -lgmp
debug_t: $(objects)
	gcc -g -Wall -Wextra -O2 source/test.c -o debug.o $(objects) -lgmp
test: $(objects)
	gcc -Wall -Wextra -O3 -static source/test.c -o test.o $(objects) -lgmp -pthread
prof: $(objects)
	gcc -pg -Wall -Wextra -O2 source/start.c -o prof.o $(objects) -lgmp -pthread
prof_t: $(objects)
	gcc -pg -Wall -Wextra -O2 source/test.c -o prof_t.o $(objects) -lgmp -pthread

prof_s: $(objects)
	gcc -pg -Wall -Wextra -O2 -static source/start.c -o prof.o $(objects) -lgmp -pthread