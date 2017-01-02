all:
	g++ -w main.cpp randomnetwork.cpp importparam.cpp connectivity.cpp cgmethod.cpp nrutil.cpp -o mikado

#SRC = $(wildcard *.cpp)
#OBJ = $(subst .cpp,.o,$(SRC))

#.PHONY = mikado test

#all: mikado

#test:
#	@echo $(SRC)
#	@echo $(OBJ)

#mikado: $(OBJ)
#	g++ $^ -o $@

#%.o: %.cpp
#	g++ -c $^ $@
