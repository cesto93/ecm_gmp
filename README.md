# How to compile and use

### Prerequisites 
- make
- gcc compiler
- gmp (The GNU Multiple Precision Arithmetic Library) [https://gmplib.org]

### Install steps
- make install
- ./start.o N B1 [B2] [MAX_ITER]
	- B2: if omitted B2 = B1 * 100
	- MAX_ITER: if omitted MAX_ITER = 2000

### Special configuration 
Edit file m_ellcurv_fact.h, you should find:
```
#define NO_LINUX => disable linux only function	
#define NO_THREAD => disable multithread 
#define THREAD_NUM => the number of threads (must be > 0) (with NO_THREAD its ignored)
```
