# RANDOMIZED SINGLE VALUE DECOMPOSITION

Contributors:

- Cenzato Matteo
- Pisante Giuseppe
- Procaccio Arianna

## GENERAL INFORMATION

This repository consists in a collection of programs that calulate some decompositions of (in general overdetermined) matrices.

Namely it can compute the [QR](https://en.wikipedia.org/wiki/QR_decomposition) and [SVD](https://en.wikipedia.org/wiki/Singular_value_decomposition) decompositions using various techniques, even randomized ones.

## MATRIX FORMAT

In order to provide an input matrix for analysis the format should be a *.txt* file where the first line contains the dimension of the matrix (**rows x cols**) and the rest are the entries of a matrix row by row with a ' ' between every column.

## COMPILING THE TESTS

In order to compile the tests we can execute the command `make all` from the root directory of the project.

This will also automatically create a folder `./build` if it is not already present where all the executables are stored.

The test files can also be generated indipendently, in fact the makefile supports also the following commands:

+ `make fullMatrix`
+ `make svd`
+ `make qr`
+ `make help`      (if for some reason you forget about them)

Also each command supports the additional options `parallel=on` and `eigen=on` which respectively compile the program using OpenMp and the matrices provided by the [Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) library (make sure to have it installed!!).

## RUNNING THE TESTS

# FullMatrix

The `fullMatrix` tests computes the most expensive operation in the whole project: the **matrix-matrix multiplication** using [Hilbert Matrices](https://en.wikipedia.org/wiki/Hilbert_matrix). In particular it takes two matrices $` A \in \mathbb{R}^{m \times n} , B \in \mathbb{R}^{n \times q} `$ and it outputs the timing results of $` C=AB , C \in \mathbb{R}^{m \times q} `$ .

The command to run it is:

```
./fullMatrix m n q
```
In particular based on the options specified at compile time (specified above) it can provide additional information and comparisons.

# SVD

The `SVD` tests computes the Singular Value Decomposition of a matrix both using the Power Method and the rSVD algorithm. In particular, it takes the matrix $` A \in \mathbb{R}^{m \times n} `$ (which can be provided by the user or is a predefined one - matrix2.txt in the test_matrices folder) and outputs the error and the timing results using two different algorithms based on the power method and one on the rSVD algorithm.

The command to run it is:
```
./svd filename.txt
```
filename.txt is optional.
The test returns also the computed matrices of the algorithm based on the power method (only of the one that performs better) and of the rSVD.

# QR

The `QR` tests computes the QR factorization both with Givens rotations and Householder algorithm. In particular it takes the matrix $` A \in \mathbb{R}^{m \times n} `$ and it outputs the timing results of the two algorithms both in serial and in parallel, showing the speed up achieved.

The command to run it is:

```
./qr m n
```

Such test also returns 4 files with the matrixes $` R \in \mathbb{R}^{m \times n} , Q \in \mathbb{R}^{m \times m} `$ computed by the 2 algorithms both serially and in parallel. Dimension on the test matrix can be added too.
Files names are printed on the output.
