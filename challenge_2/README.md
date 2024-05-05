# Compilation
## Debug mode
I set up a way to print some useful message only if you are in a debugging-phase. To use it, add the following line in the code of interest:
``` C++
DEBUG_MSG("Message to print")
```
and compile using:
```bash
make debug
```

In this way, I found useful to switch quickly from an experimental version of the code to an optimized one.

## Time testing mode
To test the time of the program use:
```bash
make optimized
```
By this way the program use the optimized version of the code.

# Usage
After compiling the program, the call is:
```bash
./main matrix_file.mtx
```

# Performance - Matrix-Vector multiplication
I tested the performance of the program using this [matrix](https://math.nist.gov/MatrixMarket/data/Harwell-Boeing/lns/lnsp_131.html). It is a matrix with 131 rows and 131 columns.
I repeated the test 10 times using the following command: 
```Bash
seq 10 | xargs -Iz ./main test/lnsp_131.mtx  
```
All test were executed on a laptop with the following characteristics:
- Model: 82UT Yoga Slim 7 Pro 14IAH7
- OS: Arch Linux x86_64
- CPU: 12th Gen Intel i7-12700H
- RAM: 16GB

I took the average of the time (in **microseconds**) spent by the program to execute the matrix-vector moltiplication:

| Test | Compressed | Decompressed |
|:----:|:--------:|:----------:|
| 1    |  1 |  299  |
| 2    |  2 |  676  |
| 3    |  2 |  646  |
| 4    |  1 |  447  |
| 5    |  1 |  589 |
| 6    |  1 |  296  |
| 7    |  1 |  299  |
| 8    |  1 |  311  |
| 9    |  2 |  698 |
| 10   |  1 |  438  |

The average time spent in matrix-vector multiplication is:
|matrix mode| mean time (microseconds)|
|:----:|:----:|
| compressed | 1,3 |
| uncompressed| 469,9 |

# Side notes
### TODOs
It may happen that there are comments in the code that begin with “TODO:”. It is a reminder for me to implement some features or to fix some bugs.

### Test_matrix.mtx
Very small matrix to see if the function *read_matrix_MM* and *operator** read correctly the matrix. 

### Terms of multiplication
To speed-up the development and debugging phase, the vector and the matrix for which I multiply are plenty of 1s. 

# Further Improvements
- Add the possibility to choose runtime if save as columnMajor or rowMajor the matrix; (It probably decrease perfomance)
- Add the possibility to pass as parameter to main the vector\matrix to multiply;
- save the resulting of that multiplication in Matrix Market Format 
