# Dependencies
This project use the following project so, before use this program, it is necessary to install it.

## muParser
To use this program you need to have installed muParser. This is the [guide](https://github.com/beltoforion/muparser/blob/master/Install.txt) on how to install it.

## MPI and OpenMP
Guides on how to download them are easily found online.

# How to use it

## Set parameters
You can set two parameters:
- `n_max`: the maximum number of iterations to reach convergence;
- `tolerance`: the tolerance to reach convergence.
- `u`: right solution to test the code; at the moment the code is set to work with the test function `u(x,y) = sin(pi*x)*sin(pi*y)` which is the solution of ```-Laplace(u) = 2*pi*pi*sin(pi*x)*sin(pi*y)```.

## Building
Make commands are: 
- debug: compile the code with flags which helps to debug; 
- mpi: compile the code with flags which optimize it;
- test_compile: compile the code in an optimized way and enable test features (**Remember** to modify the exact solution into parameters.hpp file).;
- clean: to clean all object files and executables;

## Running
To run the program you need to use the following command:
```bash
mpirun -np [number of processes] ./main [dimension of the matrix] [function f to work with] [number of openMP threads] [boundaries function]
```

As output on the screen will be printed: 
- the number of iterations to reach convergence; 
- time taken to all updates of the mesh; 
- mean time of each single update;
- [Only in Test Mode] the error between the solution and the exact solution;

While on the folder vtk_files will be saved the vtk files to visualize the mesh in a vtk reader like Paraview.

An example, similar to commands inside the test.sh script, is:
```bash
mpirun -np 4 ./main 4 "2*pi*pi*sin(pi*x)*sin(pi*y)" 4 "0"
```

## Testing
To test the program you can use the following command:
```bash
bash test.sh
```
It will ask 2 number: 
- the first one is the maximum power of 2 to use as dimension of the matrix;
- the second one is the number of processes to use (Max number of processes of your machine that MPI can use, usually is the number of cores);

The script will firstly sequentially solve the problem from a dimension of 4x4 to 2^max_power x 2^max_power and then will solve the problem in parallel using each possible configuration of processes from 2 to input of the program for each possible dimension (from 4x4 to 2^max_power x 2^max_power). Since the program test multiple configurations, output will be saved in a file called `report.txt`.

# Results
My computations time are stored in the file `report.txt`. While the file `hw.info` contains the information about the machine I used to run the tests.

# Code Organization
I have created 3 class to better divide the work and each one has its own duties to better organize the code and keep it maintainable.
The classes are:
1) `mesh_data_class`: to store all the data related to the mesh and its update; 
2) `Mesh`: to operate on the mesh like updating it, calculating error between two iterations and so on; 
3) `Solver`: to solve the problem using Jacobi method leaning on previous classes.

In this way possible future improvements and/or changes will be easier to implement.

# TODO
- Other methods to solve the problem, like Schwarz iteration;
- Adding more tricky way to insert boundaries conditions, like Newmann or Robin boundaries;
- Reading the input of the program from a json file;