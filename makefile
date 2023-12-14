#My flags
OPTIMIZATION=-O2 -march=native
WARNINGS= -Wall -pedantic
EIGEN=-I${mkEigenInc}
SRC=./src
BUILD=./build
FULLMATRIX_MAIN=fullMatrixTest.cpp
FULLMATRIX_OBJS=$(FULLMATRIX_MAIN:.cpp=.o)
FULLMATRIX_OUTS=$(FULLMATRIX_MAIN:.cpp=)

ifdef parallel
CXXFLAGS+=-fopenmp
LDFLAGS+=-fopenmp
else
CXXFLAGS+=-Wno-unknown-pragmas
endif

#Compiler version
CXX=g++
#Preprocessor flags
CPPFLAGS+= -I./include
#Compiler flags
CXXFLAGS+= -std=c++20 $(WARNINGS)
#Linker flags
LDFLAGS+= $(OPTIMIZATION)

#Make directives that have no dependencies
.PHONY= all optimized help clean distclean

help:
	@echo "All compilation commands: "
	@echo "  all         (default)"
	@echo "  fullMatrix           "
	@echo "  clean                "
	@echo "  distclean             "

all:
	g++  src/main.cpp -o build/main $(CXXFLAGS) $(WARNINGS) $(CPPFLAGS)
	$(MAKE) clean

optimized:
	$(MAKE) CXXFLAGS="$(OPTIMIZATION)"

eigen:
	$(MAKE) CPPFLAGS="$(EIGEN)"

clean: 
	$(RM) $(BUILD)/*.o

distclean:
	$(RM) $(BUILD)/*

#I have three sections -> fullMatrix, svd, qr

fullMatrix: fullMatrix.o
	$(CXX) $(BUILD)/$(FULLMATRIX_OBJS) -o $(BUILD)/$(FULLMATRIX_OUTS) $(LDFLAGS)
	$(MAKE) clean

fullMatrix.o: $(SRC)/$(FULLMATRIX_MAIN)
	$(CXX) $(SRC)/$(FULLMATRIX_MAIN) -c -o $(BUILD)/$(FULLMATRIX_OBJS) $(CXXFLAGS)
