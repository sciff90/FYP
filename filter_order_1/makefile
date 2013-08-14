#gaussian Makefile
#Define Objects

OBJ = obj
SRC = src
PLOTS = plots
DATA = data
EXEC = xfilter_order_1

OBJECTS = ${OBJ}/filter_order_1.o

#Define Compiler

#############gcc############
COMP = gcc
	
	#Define Switches
	#SWITCH = -O2
	SWITCH = -g


LIB_DIR = 
LIB = -lgsl -lgslcblas -lm
IDIR = -I.
##############################


#Makefile
all: $(OBJECTS)
		$(COMP) $(SWITCH) ${OBJECTS} $(IDIR) -o ${EXEC} $(LIB_DIR) $(LIB)
			
${OBJ}/filter_order_1.o: ${SRC}/filter_order_1.c
		$(COMP) -c $(SWITCH) $(LIB_DIR) $(LIB) $(IDIR) ${SRC}/filter_order_1.c -o ${OBJ}/filter_order_1.o
			
clean:
		rm ${OBJ}/*.o
		rm ${PLOTS}/*.ps
		rm ${DATA}/*.dat

#Cleaning Everything
cleanest: clean
		rm ${EXEC}
# End of the makefile
