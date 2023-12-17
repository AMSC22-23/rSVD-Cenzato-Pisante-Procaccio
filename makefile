#My flags
OPTIMIZATION=-O2 -march=native
WARNINGS= -Wall -pedantic 
EIGEN=${mkEigenInc}
BUILD=./build

#Files for the fullMatrix test
FULLMATRIX_TEST=fullMatrixTestTiming.cpp

#Files for the QR_Decomposition test
QR_TEST=QR_DecompositionTest.cpp QR_Decomposition.cpp QR_Decomposition_parallel.cpp

#Files for the SVD test
SVD_TEST=svd_test.cpp svd.cpp QR_Decomposition.cpp

ifdef parallel
CXXFLAGS+=-fopenmp
LDFLAGS+=-fopenmp
else
CXXFLAGS+=-Wno-unknown-pragmas
endif

ifdef eigen
CXXFLAGS+=-I$(EIGEN)
LDFLAGS+=-I$(EIGEN)
CPPFLAGS+=-DEIGEN
endif

#Compiler version
CXX=g++
#Preprocessor flags
CPPFLAGS+= -I./include 
#Compiler flags
CXXFLAGS+= -std=c++20 $(WARNINGS)
#Linker flags
LDFLAGS+= $(OPTIMIZATION)
#Path for the files 
VPATH=./src

#Make directives that have no dependencies
.PHONY= all optimized help clean distclean

help:
	@echo "All compilation commands: "
	@echo "  all         (default)"
	@echo "  fullMatrix           "
	@echo "  clean                "
	@echo "  distclean            "

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

fullMatrix: $(FULLMATRIX_TEST)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS)
	$(MAKE) clean

qr: $(QR_TEST)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) 
	$(MAKE) clean

svd: $(SVD_TEST)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) 
	$(MAKE) clean

#$(QR_TEST):  
#	$(CXX) $@ -c -o $< $(CXXFLAGS) $(CPPFLAGS) -I${mkEigenInc}

#%.o: %.cpp
#	$(CXX) $< -c -o $@ $(CXXFLAGS) $(CPPFLAGS)
