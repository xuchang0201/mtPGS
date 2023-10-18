# -----------------------------------------------------------------
#   Makefile for mtPGS
# ---------------------------------------------------------------------

# Set the file type
OUTPUT1 = mtPGS_int_ext
OUTPUT2 = mtPGS_int_only

# Put C++ complier 
CXX = g++

# Set complier flags 
CXXFLAG = -O2 -std=c++11 -lm -Wall -larmadillo -fopenmp 
all: $(OUTPUT1) $(OUTPUT2)
$(OUTPUT1): main_software_int_ext.o mtPGS_function.o mtPGS_model_int_ext.o mtPGS_int_ext.o
	$(CXX) main_software_int_ext.o mtPGS_function.o mtPGS_model_int_ext.o mtPGS_int_ext.o -o $(OUTPUT1) $(CXXFLAG) 
$(OUTPUT2): main_software_int_only.o mtPGS_function.o mtPGS_model_int_only.o mtPGS_int_only.o 
	$(CXX) main_software_int_only.o mtPGS_function.o mtPGS_model_int_only.o mtPGS_int_only.o -o $(OUTPUT2) $(CXXFLAG) 
main_software_int_ext.o: main_software_int_ext.cpp
	$(CXX) -c main_software_int_ext.cpp
main_software_int_only.o: main_software_int_only.cpp
	$(CXX) -c main_software_int_only.cpp 
mtPGS_model_int_ext.o: mtPGS_model_int_ext.cpp mtPGS_model_int_ext.hpp
	$(CXX) -c mtPGS_model_int_ext.cpp $(CXXFLAG)
mtPGS_model_int_only.o: mtPGS_model_int_only.cpp mtPGS_model_int_only.hpp
	$(CXX) -c mtPGS_model_int_only.cpp $(CXXFLAG)
mtPGS_function.o: mtPGS_function.cpp mtPGS_function.hpp 
	$(CXX) -c mtPGS_function.cpp $(CXXFLAG)
mtPGS_int_ext.o: mtPGS_int_ext.cpp mtPGS_int_ext.hpp
	$(CXX) -c mtPGS_int_ext.cpp $(CXXFLAG)
mtPGS_int_only.o: mtPGS_int_only.cpp mtPGS_int_only.hpp
	$(CXX) -c mtPGS_int_only.cpp $(CXXFLAG)

clean:
	rm -f *.o  mtPGS_int_ext mtPGS_int_only
