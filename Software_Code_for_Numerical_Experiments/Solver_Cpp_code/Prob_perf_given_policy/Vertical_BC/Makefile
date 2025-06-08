CC = g++
CXXFLAGS = -std=c++17 -O3 -fopenmp
BOOSTPATH = ../../cpp/boost_1_79_0
EIGENPATH = ../../cpp/eigen-3.4.0/
INCLUDE = -I $(BOOSTPATH) -I $(EIGENPATH) -I ./include/
LIB = -L $(BOOSTPATH)
SRC = ./src/
TEST = 

CPPFILES = $(SRC)main.cpp \
					$(SRC)CGrid_1D.cpp \
					$(SRC)CGrid_2D.cpp \
					$(SRC)CPolicyEvaluation.cpp \
					$(SRC)CSolver.cpp \
					
main: $(CPPFILES)
	$(CC) $(CXXFLAGS) -o $@ $^ $(INCLUDE) $(LIB)

run: main
		./main $(TEST)

debug: $(CPPFILES)
	g++ -o1 -g -std=c++17 $(INCLUDE) $^ $(LIB)
	gdb -ex=r --args 1 $(TEST)

clean:
	rm *.dat main
