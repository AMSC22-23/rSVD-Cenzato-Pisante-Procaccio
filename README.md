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

TODO

# QR

TODO
