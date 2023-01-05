# How to compile and use

### Prerequisites 
- make
- gcc compiler
- gmp (The GNU Multiple Precision Arithmetic Library) [https://gmplib.org]

### Install steps
- make
- ./start.o N B1 [B2] [MAX_ITER]
	- B2: if omitted B2 = B1 * 100
	- MAX_ITER: if omitted MAX_ITER = 2000

### Special configuration 
- make mm: use mm_ellcurv lib
```
#define THREAD_NUM => the number of threads (must be > 0) (with NO_THREAD its ignored)
```
