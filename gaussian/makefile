#gaussian Makefile
#Define Objects

OBJ = obj
SRC = src
PLOTS = plots
DATA = data
EXEC = xgaussian

OBJECTS = ${OBJ}/gaussian.o

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
all: DIR $(OBJECTS)
		$(COMP) $(SWITCH) ${OBJECTS} $(IDIR) -o ${EXEC} $(LIB_DIR) $(LIB)
			
${OBJ}/gaussian.o: ${SRC}/gaussian.c
		$(COMP) -c $(SWITCH) $(LIB_DIR) $(LIB) $(IDIR) ${SRC}/gaussian.c -o ${OBJ}/gaussian.o
DIR:
		mkdir -p $(OBJ) $(PLOTS) $(DATA)
			
clean:
		rm -rf ${OBJ}
		rm -rf ${PLOTS}
		rm -rf ${DATA}

#Cleaning Everything
cleanest: clean
		rm -f ${EXEC}
# End of the makefile
