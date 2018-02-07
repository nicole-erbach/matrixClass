# set compiler
CC = g++

# compiler flags: 
CFLAGS = -Wall -g -std=c++11 -fopt-info -I /usr/include/eigen3 -I incl/

# List with source code files
SRC = matrix.cpp test.cpp

# Libs to link
LIB = -lblas

# For list with object files, replace .cpp ending with .o ending
OBJ = $(SRC:%.cpp=%.o)

# main target: build all .so files to build/
main: $(OBJ)
	  $(CC) $(OBJ) $(LIB) -o $(OUTPUT)

OUTPUT = test.out

# Generic rule to create .o files from .cpp files
%.o: %.cpp
	$(CC) $(CFLAGS) -c $< -o $@

# target to print variables: make print-VARIABLE
# print-%  : ; @echo $* = $($*)

.PHONY: clean

# Target that removes all .so files
clean:
	rm -rf $(OBJ) $(OUTPUT)
