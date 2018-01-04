INCLUDE = -I/usr/include/
LIBDIR  = -L/usr/lib

FLAGS = -Wall
CC = g++                                  # change to gcc if using C
CFLAGS = $(FLAGS) $(INCLUDE)
LIBS =  -lglut -lGL -lGLU -lGLEW -lm

#All: clean eigenSolver.o visual.o EigenProj                             # change your_app.

EigenProj: main.o visual.o 
	$(CC) $(CFLAGS) -o $@ $(LIBDIR) $^ $(LIBS) # The initial white space is a tab

main.o: main.cpp visual.h
	$(CC) -c main.cpp
visual.o: visual.cpp visual.h
	$(CC) -c visual.cpp



clean:  
	rm -rf *.o  
	rm -rf EigenProj  
