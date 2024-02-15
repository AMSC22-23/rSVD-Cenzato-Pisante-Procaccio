#My flags
OPTIMIZATION=-O3 -march=native
WARNINGS= -Wall -pedantic 
EIGEN=${mkEigenInc}
BUILD=./build

#Files for the fullMatrix test
FULLMATRIX_TEST=fullMatrixTestTiming.cpp

#Files for the QR_Decomposition test
QR_TEST=QR_DecompositionTest.cpp QR_Decomposition.cpp QR_Decomposition_parallel.cpp

#Files for the SVD test
SVD_TEST=svd_test.cpp svd.cpp QR_Decomposition_parallel.cpp

#Files for image compression
IMAGE_TEST=svd.cpp QR_Decomposition_parallel.cpp imagecompr.cpp

#Files for pca 
PCA_TEST=cancer.cpp svd.cpp QR_Decomposition_parallel.cpp

#Files for the benchmarks
BENCHMARKS=benchmarkTimings.cpp

ifdef parallel
CXXFLAGS+=-fopenmp
#LDFLAGS+=-fopenmp
else
CXXFLAGS+=-Wno-unknown-pragmas
endif

ifdef eigen
CXXFLAGS+=-I$(EIGEN)
LDFLAGS+=-I$(EIGEN)
CPPFLAGS+=-DEIGEN
endif

ifdef RGB
CPPFLAGS+=-DRGB
endif

ifdef lazy
CPPFLAGS+=-DLAZY
endif

#Compiler version
CXX=mpic++
#Preprocessor flags
CPPFLAGS+= -I./include 
#Compiler flags
CXXFLAGS+= -std=c++20 $(WARNINGS) $(OPTIMIZATION)
#Linker flags
#LDFLAGS+= $(OPTIMIZATION)
#Path for the files 
VPATH=./src

#Make directives that have no dependencies
.PHONY= all help clean distclean

help:
	@echo "All compilation commands: "
	@echo "  all         (default)"
	@echo "  fullMatrix           "
	@echo "  qr                   "
	@echo "  svd                  "
	@echo "  compression		      "
	@echo "  pca				          "
	@echo "  clean                "
	@echo "  distclean            "
	@echo "  benchmarks           "

all:
	mkdir -p $(BUILD)
	$(MAKE) qr
	$(MAKE) svd
	$(MAKE) fullMatrix
	$(MAKE) compression
	$(MAKE) pca
	$(MAKE) benchmarks

clean: 
	$(RM) $(BUILD)/*.o

distclean:
	$(RM) $(BUILD)/*

#I have three sections -> fullMatrix, svd, qr

fullMatrix: $(FULLMATRIX_TEST)
	mkdir -p $(BUILD)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) -DBYROWS
	$(MAKE) clean

qr: $(QR_TEST)
	mkdir -p $(BUILD)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) 
	$(MAKE) clean

svd: $(SVD_TEST)
	mkdir -p $(BUILD)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) 
	$(MAKE) clean

compression: $(IMAGE_TEST)
	mkdir -p $(BUILD)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) 
	$(MAKE) clean

pca: $(PCA_TEST)
	mkdir -p $(BUILD)
	$(CXX) $^ -o $(BUILD)/$@ $(LDFLAGS) $(CPPFLAGS) $(CXXFLAGS) 
	$(MAKE) clean

benchmarks: $(BENCHMARKS)
	mkdir -p $(BUILD)
	$(CXX) $^ -o $(BUILD)/b_lazyprodlong  $(LDLFLAGS) $(CPPFLAGS) $(CXXFLAGS) -DLAZYPRODLONG
	$(CXX) $^ -o $(BUILD)/b_lazyprodshort $(LDLFLAGS) $(CPPFLAGS) $(CXXFLAGS) -DLAZYPRODSHORT
	$(CXX) $^ -o $(BUILD)/b_lazysumlong   $(LDLFLAGS) $(CPPFLAGS) $(CXXFLAGS) -DLAZYSUMLONG
	$(CXX) $^ -o $(BUILD)/b_lazysumshort  $(LDLFLAGS) $(CPPFLAGS) $(CXXFLAGS) -DLAZYSUMSHORT
	$(CXX) $^ -o $(BUILD)/b_lazymixed     $(LDLFLAGS) $(CPPFLAGS) $(CXXFLAGS) -DMIXED
	$(CXX) $^ -o $(BUILD)/b_matmult       $(LDLFLAGS) $(CPPFLAGS) $(CXXFLAGS) 
	$(MAKE) clean

#$(QR_TEST):  
#	$(CXX) $@ -c -o $< $(CXXFLAGS) $(CPPFLAGS) -I${mkEigenInc}

#%.o: %.cpp
#	$(CXX) $< -c -o $@ $(CXXFLAGS) $(CPPFLAGS)
