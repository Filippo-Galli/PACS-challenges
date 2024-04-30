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
| 1    |  20 |  8999  |
| 2    |  22 |  14126  |
| 3    |  23 |  9653  |
| 4    |  21 |  10522  |
| 5    |  17 |  8158 |
| 6    |  19 |  9739  |
| 7    |  26 |  8287  |
| 8    |  17 |  8151  |
| 9    |  22 |  8689 |
| 10   |  20 |  8805  |

The average time spent in matrix-vector multiplication is:
|matrix mode| mean time (microseconds)|
|:----:|:----:|
| compressed | 20,7 |
| uncompressed| 9512 |

## Without dimension check
Without the dimension check the average time spent in matrix-vector multiplication is:

| Test | Compressed | Decompressed |
|:----:|:--------:|:----------:|
| 1    |  14 |  5684  |
| 2    |  12 |  4523  |
| 3    |  13 |  4744  |
| 4    |  11 |  3912  |
| 5    |  11 |  3798 |
| 6    |  11 |  3913  |
| 7    |  15 |  4220  |
| 8    |  13 |  3954  |
| 9    |  15 |  5174 |
| 10   |  9 |  3960  |

The average time spent in matrix-vector multiplication is:
|matrix mode| mean time (microseconds)| % of improvement |
|:----:|:----:|:----:|
| compressed | 12,4 | -40,1% |
| uncompressed| 4388 | -53,87% |

# Side notes
### TODOs
It may happen that there are comments in the code that begin with “TODO:”. It is a reminder for me to implement some features or to fix some bugs.

### Test_matrix.mtx
Very small matrix to see if the function *read_matrix_MM* read correctly the matrix. 


