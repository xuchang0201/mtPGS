# -----------------------------------------------------------------
#   Makefile for mtPGS
# ---------------------------------------------------------------------

# Set the file type
OUTPUT = mtPGS

# Put C++ complier 
CXX = g++

# Set complier flags 
CXXFLAG = -O2 -std=c++11 -lm -Wall -larmadillo -fopenmp 
all: $(OUTPUT) 
$(OUTPUT): main_software.o mtPGS_function.o mtPGS_model.o mtPGS.o
	$(CXX) main_software.o mtPGS_function.o mtPGS_model.o mtPGS.o -o $(OUTPUT) $(CXXFLAG) 
main_software.o: main_software.cpp
	$(CXX) -c main_software.cpp 
mtPGS_model.o: mtPGS_model.cpp mtPGS_model.hpp
	$(CXX) -c mtPGS_model.cpp $(CXXFLAG)
mtPGS_function.o: mtPGS_function.cpp mtPGS_function.hpp 
	$(CXX) -c mtPGS_function.cpp $(CXXFLAG)
mtPGS.o: mtPGS.cpp mtPGS.hpp
	$(CXX) -c mtPGS.cpp $(CXXFLAG)

clean:
	rm -f *.o  mtPGS 
