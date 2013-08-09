#mcmc_serial Makefile
#Define Objects

OBJ = obj
SRC = src
PLOTS = plots
DATA = data
EXEC = xmcmc_serial

OBJECTS = ${OBJ}/mcmc_serial.o

#Define Compiler

#############gcc############
COMP = gcc
	
	#Define Switches
	#SWITCH = -O2
	SWITCH = -g


LIB_DIR = 
LIB = -lgsl -lgslcblas
IDIR = -I.
##############################


#Makefile
all: $(OBJECTS)
		$(COMP) $(SWITCH) ${OBJECTS} $(IDIR) -o ${EXEC} $(LIB_DIR) $(LIB)
			
${OBJ}/mcmc_serial.o: ${SRC}/mcmc_serial.c
		$(COMP) -c $(SWITCH) $(LIB_DIR) $(LIB) $(IDIR) ${SRC}/mcmc_serial.c -o ${OBJ}/mcmc_serial.o
			
clean:
		rm ${OBJ}/*.o
		rm ${PLOTS}/*.ps
		rm ${DATA}/*.dat

#Cleaning Everything
cleanest: clean
		rm ${EXEC}
# End of the makefile
