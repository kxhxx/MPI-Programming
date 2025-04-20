# MPI-Programming
Basic codes in C++ for parallel computing using MPI.


1. Monte Carlo Method for Estimating Pi
   
Description:
This program estimates π using the Monte Carlo method, where points are randomly generated, and the ratio of points inside a unit circle is used to estimate π.

Requirements:
MPI library (OpenMPI or MPICH),  
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ monte_carlo_pi.cpp -o monte_carlo_pi

mpirun -np <num_processes> ./monte_carlo_pi


2. Heat Distribution Simulation
   
Description: 
This MPI program simulates heat distribution across a grid using a numerical solver.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_heat_simulation.cpp -o mpi_heat_simulation

mpirun -np <num_processes> ./mpi_heat_simulation


3. Matrix Multiplication
   
Description:
This MPI program multiplies two 70×70 matrices in parallel. It also measures execution time and compares it with a sequential implementation.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_matrix_multiply.cpp -o mpi_matrix_multiply

mpirun -np <num_processes> ./mpi_matrix_multiply


4. Parallel Dot Product
   
Description:
This program computes the dot product of two large vectors in parallel using MPI.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_dot_product.cpp -o mpi_dot_product

mpirun -np <num_processes> ./mpi_dot_product


5. Parallel Matrix Transposition
   
Description:
This MPI program transposes a matrix in parallel using row-wise and column-wise communication.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_matrix_transpose.cpp -o mpi_matrix_transpose

mpirun -np <num_processes> ./mpi_matrix_transpose


6. Parallel Sorting (Odd-Even Sort)
   
Description:
This program sorts an array using Odd-Even Sort in parallel using MPI.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_odd_even_sort.cpp -o mpi_odd_even_sort

mpirun -np <num_processes> ./mpi_odd_even_sort

 
7. Parallel Prefix Sum (Scan)
   
Description: 
This MPI program performs a parallel prefix sum (scan) operation.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_prefix_sum.cpp -o mpi_prefix_sum

mpirun -np <num_processes> ./mpi_prefix_sum


8. Parallel Reduction
   
Description:
This MPI program demonstrates parallel reduction, where multiple values are combined into a single result using MPI_Reduce.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_reduction.cpp -o mpi_reduction

mpirun -np <num_processes> ./mpi_reduction

9. Estimating Pi using MPI Bcast and MPI Reduce
    
Description:
This MPI program approximates the value of π using numerical integration.
It follows these steps:
The number of steps (num_steps) is broadcasted using MPI_Bcast.
Each process computes its partial sum.
The partial sums are combined using MPI_Reduce.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_pi.cpp -o mpi_pi

mpirun -np <num_processes> ./mpi_pi

10. DAXPY Loop
    
Description:
This MPI program performs the DAXPY operation, where each element of the vector X is updated as:

𝑋[𝑖]=𝑎×𝑋[𝑖]+𝑌[𝑖]

It also measures the speedup of MPI execution compared to a sequential implementation.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ daxpy_mpi.cpp -o daxpy

mpirun -np <num_processes> ./daxpy


11. Finding Prime Numbers
    
Description:
This MPI program finds prime numbers up to a given maximum using a master-slave model.
The master process assigns numbers to test.
The slave processes test numbers and report results back.

Requirements:
MPI library (OpenMPI or MPICH), 
C++ compiler supporting MPI (mpic++)

Compilation & Execution:

mpic++ mpi_primes.cpp -o mpi_primes

mpirun -np <num_processes> ./mpi_primes <max_value>


General Notes:

Ensure MPI is installed (mpicc --version).

Use mpirun -np <num_processes> to run the programs.

Modify matrix size, number of iterations, or other parameters in the source code for different test cases.


12. CUDA Threads Performing Different Tasks
    
Description:
This program demonstrates the use of CUDA threads to perform distinct tasks in parallel. Thread 0 calculates the sum of the first N integers using an iterative loop, while Thread 1 calculates the same sum using the direct mathematical formula.

Requirements:
NVIDIA GPU with CUDA support
CUDA Toolkit
C++ compiler with nvcc (NVIDIA CUDA compiler)

Compilation & Execution:

nvcc different_tasks.cu -o different_tasks

./different_tasks


13. Parallel Merge Sort (CPU Pipelined)
    
Description:
This program implements a parallel merge sort using standard C++ threads. The array is recursively divided and sorted in parallel before being merged in a pipelined fashion.

Requirements:
C++11 or newer
Thread support (-pthread)
Any standard C++ compiler (e.g., g++, clang++)

Compilation & Execution:

g++ -std=c++11 -O3 merge_sort_pipelined.cpp -o merge_sort_cpu -pthread

./merge_sort_cpu


14. Parallel Merge Sort (CUDA)
    
Description:
This CUDA-based implementation of a parallel merge sort uses a simplified bitonic sorting approach to perform sorting on the GPU using many threads.

Requirements:
NVIDIA GPU with CUDA support
CUDA Toolkit
C++ compiler with nvcc

Compilation & Execution:

nvcc -O3 merge_sort_cuda.cu -o merge_sort_cuda

./merge_sort_cuda


15. CUDA Vector Addition with Profiling:
    
Description:
This CUDA program adds two vectors of size N using statically defined device memory (no cudaMalloc). It includes:
Kernel timing using CUDA Events
Device query to compute theoretical memory bandwidth
Runtime-based measured bandwidth calculation

Requirements:
NVIDIA GPU with CUDA support
CUDA Toolkit (10.0+)
nvcc compiler

Optional: nvprof or nsys for profiling

Compilation & Execution: 
nvcc -O3 vector_add_static.cu -o vector_add

./vector_add
Profiling (Optional)

nvprof ./vector_add
or
nsys profile ./vector_add
