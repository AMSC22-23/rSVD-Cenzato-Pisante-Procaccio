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
+ `make pca`
+ `make compression`
+ `make benchmarks`
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

The `SVD` test computes the Singular Value Decomposition of a matrix both using the Power Method and the rSVD algorithm. In particular, it takes a gaussian matrix $` A \in \mathbb{R}^{m \times n} `$ (where m and n can be provided by the user) and outputs the error and the timing results of the different algorithms.

The command to run it is:
```
./svd [flag] [m] [n] [r]
```
Where flag, m, n and r are optional:

+ `flag` = 0 (which is the default): computes 2 algorithms that use the power method and the rSVD; `flag` = 1: computes PM1 and PM2; `flag` = 2: computes PM1 and rSVD; `flag` = 3: computes PM1 and pseudo-inverse using PM1; any other value computes only PM1.

+ `m` and `n` are used to set the dimensions of the matrix (default are m = n = 10).

+ `r` is used to set the target rank in the rSVD (default is r = 5).


# QR

The `QR` tests computes the QR factorization both with Givens rotations and Householder algorithm. In particular it takes the matrix $` A \in \mathbb{R}^{m \times n} `$ and it outputs the timing results of the two algorithms both in serial and in parallel, showing the speed up achieved.

The command to run it is:

```
./qr m n
```

Such test also returns 4 files with the matrixes $` R \in \mathbb{R}^{m \times n} , Q \in \mathbb{R}^{m \times m} `$ computed by the 2 algorithms both serially and in parallel. Dimension on the test matrix can be added too.
Files names are printed on the output.


# PCA

In this test, we performed the principal component analysis (PCA) on a dataset from the FDA-NCI Clinical Proteomics Program Databank.

Each column of the dataset represents measurements taken from a patient. There are 216 columns
representing 216 patients, out of which 121 are ovarian cancer patients and 95 are normal patients.
Each row represents the ion intensity level at a specific mass-charge value indicated in MZ. There
are 2000 mass-charge values, and each row represents the ion-intensity levels of the patients at that
particular mass-charge value.

The command to run the test is:

```
./pca
```

It returns a matrix containing the first 50 principal components of the dataset.

# BENCHMARKS

To view the effectiveness of the lazy evaluation compared to a naive implementation, we can run the command `make benchmarks` which creates the executables containing some default operations explained in detail in the report.

In order to run any of the tests the command is :

```
./b_matmult start end step
```

Where `start` is the starting dimension of the matrix, `end` the last dimension, `step` how much the matrix increases in dimension at each iteration (it is additive, not multiplicative).

**NB**: Remember to add `make ... lazy=on` in order to use lazy evaluation.
        Also when compiling with Eigen it gives a wall of warings because of the deprecated functions used in the Eigen library when compiling with `std=c++20`, so nothing can be done about it.
