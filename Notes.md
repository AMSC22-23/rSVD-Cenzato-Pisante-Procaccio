# Some comments #

In QR_Decomposition.cpp there is a dangerous error linked to the use of `abs()` Unfortunately, due to pollotion with C headers in the current implementation of the standard library in GNU compilers, you have sometimes at disposal
` int abs(int)` that works **on integers**. If you want the general `abs()` that works on any type you must use `std::abs()` that is defined in the `cmath` header. Not `abs()` **you should correct this error**. 


In general, use the full qualified name on all standard library functions and classes. You avoid mistakes.

Readme.md file is poor! Please complete it with a description of the code and the algorithms used. You can also add a section with the compilation instructions.

In QR_Decomposition.cpp you have code duplication. You can avoid it using functions, even lambda expression. It will make the code more readable.

You should either have the Makefile create the `build` directory or add it to the repository. It is not a good idea to have the user create it (without knowing thet it is needed!).

For the rest, it looks a reasonable implementation. I have not tested it, hoeever. The use of modern c++ is good.

