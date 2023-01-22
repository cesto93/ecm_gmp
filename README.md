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

### Parameters to use
|digits	|B1	    |Curvs	|k-size	|k-bits|
|-------|-----------|-----------|-------|------|
|15	|2000       |25		|45	|2880  |
|20	|11000	    |90		|249	|15936 |
|25	|50000	    |300	|1127	|72128 |
|30	|250000	    |700	|5639	|360896|
|35	|1000000    |1800	|	|      |
|40	|3000000    |5100	|	|      |
|45	|11000000   |10600	|	|      |
|50	|43000000   |19300	|	|      |
|55	|110000000  |49000	|	|      |
|60	|260000000  |124000	|	|      |
|65	|850000000  |210000	|	|      |
|70	|2900000000 |340000	|	|      |
