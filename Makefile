FLAGS=-std=c++11 -g -Wall

ols: main.o #matrix.o
	g++ -o ols  main.o #matrix.o

main.o: main.cpp matrix.h
	g++ $(FLAGS) -c main.cpp

matrix.o: matrix.cpp matrix.h
	g++ $(FLAGS) -c matrix.cpp

clean: 
	rm -f ols *.o

rebuild: 
	make clean
	make
