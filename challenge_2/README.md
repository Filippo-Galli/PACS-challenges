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
I repeated the test 10 times and I took the average of the time (in **microseconds**) spent by the program to execute the matrix-vector moltiplication:

| Test | Compressed | Decompressed |
|:----:|:--------:|:----------:|
| 1    |  26 |  6628  |
| 2    |  20 |  4813  |
| 3    |  23 |  6843  |
| 4    |  25 |  6681  |
| 5    |  44 |  10195 |
| 6    |  30 |  8276  |
| 7    |  27 |  7490  |
| 8    |  22 |  8601  |
| 9    |  36 |  10757 |
| 10   |  22 |  8077  |

The average time spent in matrix-vector multiplication is:
|matrix mode| mean time (microseconds)|
|:----:|:----:|
| compressed | 27,5 |
| uncompressed|7836 |

All test were executed on a laptop with the following characteristics:
- Model: 82UT Yoga Slim 7 Pro 14IAH7
- OS: Arch Linux x86_64
- CPU: 12th Gen Intel i7-12700H
- RAM: 16GB

# Side notes
### TODOs
It may happen that there are comments in the code that begin with “TODO:”. It is a reminder for me to implement some features or to fix some bugs.

### Test_matrix.mtx
Very small matrix to see if the function *read_matrix_MM* read correctly the matrix. 


