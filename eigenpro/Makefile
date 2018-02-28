INCLUDE = -I/usr/include/
LIBDIR  = -L/usr/lib

FLAGS = -Wall
CC = g++                                  # change to gcc if using C
CFLAGS = $(FLAGS) $(INCLUDE)
LIBS =  -lglut -lGL -lGLU -lGLEW -lm

#All: clean eigenSolver.o visual.o EigenProj                             # change your_app.

EigenProj: main.o visual3d.o PoissonBuilder.o
	$(CC) $(CFLAGS) -o $@ $(LIBDIR) $^ $(LIBS) # The initial white space is a tab

main.o: main.cpp visual3d.h PoissonBuilder.h
	$(CC) -c main.cpp
visual3d.o: visual3d.cpp visual3d.h
	$(CC) -c visual3d.cpp
PoissonBuilder.o: PoissonBuilder.h
	$(CC) -c PoissonBuilder.cpp


clean:  
	rm -rf *.o  
	rm -rf EigenProj  
