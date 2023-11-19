all:
	g++ -I${mkEigenInc} src/main.cpp -o build/main -Wall

test:
	g++ src/fullMatrixTest.cpp -o build/fullMatrixTest -Wall
