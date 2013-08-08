#Random Walk 2D Makefile
#Define Objects

OBJ = obj
SRC = src
PLOTS = plots
DATA = data
EXEC = xrandom_walk_2d

OBJECTS = ${OBJ}/random_walk_2d.o

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
			
${OBJ}/random_walk_2d.o: ${SRC}/random_walk_2d.c
		$(COMP) -c $(SWITCH) $(LIB_DIR) $(LIB) $(IDIR) ${SRC}/random_walk_2d.c -o ${OBJ}/random_walk_2d.o
			
clean:
		rm ${OBJ}/*.o
		rm ${PLOTS}/*.ps
		rm ${DATA}/*.dat

#Cleaning Everything
cleanest: clean
		rm ${EXEC}
# End of the makefile
