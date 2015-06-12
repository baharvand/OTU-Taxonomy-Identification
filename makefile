CC=g++
CFLAGS=-c -g -Wall -O2 -I/usr/include/libxml2/
LIB=-lgsl -lgslcblas -lm -larmadillo -lmlpack

taco: main.o parseArgs.o OTUmatrix.o dictionary.o correlation.o node.o graph.o mrf.o
	$(CC) -g -Wall -O2 -o taco main.o parseArgs.o OTUmatrix.o dictionary.o correlation.o node.o \
								graph.o mrf.o $(LIB)

main.o: main.cpp parseArgs.h OTUmatrix.h mrf.h debug.h
	$(CC) $(CFLAGS) main.cpp

parseArgs.o: parseArgs.h parseArgs.cpp debug.h
	$(CC) $(CFLAGS) parseArgs.cpp

OTUmatrix.o: OTUmatrix.h OTUmatrix.cpp parseArgs.h dictionary.h correlation.h mrf.h debug.h
	$(CC) $(CFLAGS) OTUmatrix.cpp

dictionary.o: dictionary.h dictionary.cpp
	$(CC) $(CFLAGS) dictionary.cpp

correlation.o: correlation.h correlation.cpp
	$(CC) $(CFLAGS) correlation.cpp

node.o: node.h node.cpp
	$(CC) $(CFLAGS) node.cpp

graph.o: graph.h graph.cpp node.h dictionary.h
	$(CC) $(CFLAGS) graph.cpp

mrf.o: mrf.h mrf.cpp graph.h
	$(CC) $(CFLAGS) mrf.cpp

clean:
	rm -f taco *.o .*.swp *.log
