#mcmc CUDA Makefile
#Define Objects
OBJ = obj
SRC = src
DATA = data
EXEC = mcmc_cuda
MATLAB_DIR = /usr/local/MATLAB/R2013a

OBJECTS = $(MKDIR) ${OBJ}/mcmc_cuda.o

#Define Compiler

#############nvcc############
COMP = nvcc
	
#Define Switches
#SWITCH = -O2 -fopenmp -std=c++0x
SWITCH = -arch=sm_20 -m 64 -Xcompiler -fPIC 

#Define Libraries
#LIB_DIR = 
#LIB = 
IDIR = -I. -I $(MATLAB_DIR)/extern/include -I $(MATLAB_DIR)/toolbox/distcomp/gpu/extern/include
##############################


#Makefile
all: DIR $(OBJECTS)	
		mex  -cxx ${OBJECTS} 
			
${OBJ}/mcmc_cuda.o: ${SRC}/mcmc_cuda.cu
		$(COMP) -c $(SWITCH) $(LIB_DIR) $(LIB) $(IDIR) ${SRC}/mcmc_cuda.cu -o ${OBJ}/mcmc_cuda.o
DIR:
		mkdir -p $(OBJ) $(DATA)	
clean:
		rm -rf ${OBJ}		
		rm -rf ${DATA}

#Cleaning Everything
cleanest: clean
		rm -f ${EXEC}.mexa64
# End of the makefile
