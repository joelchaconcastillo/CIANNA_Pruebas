main: MA.o MPP.o Model.o localsearches.o MPP_Problem.o recomb.o global.o main.cpp 
	g++ -std=c++11 -O2 -o $@ main.cpp MA.o MPP.o Model.o localsearches.o MPP_Problem.o recomb.o global.o

MPP.o: MPP.h MPP.cpp
	g++ -std=c++11 -O2 -c -o $@ MPP.cpp

Model.o: MPP.h Model.cpp
	g++ -std=c++11 -O2 -c -o $@ Model.cpp

recomb.o: MPP.h recomb.cpp
	g++ -std=c++11 -O2 -c -o $@ recomb.cpp

MPP_Problem.o: MPP_Problem.h MPP_Problem.cpp
	g++ -std=c++11 -O2 -c -o $@ MPP_Problem.cpp

localsearches.o: MPP.h localsearches.cpp
	g++ -std=c++11 -O2 -c -o $@ localsearches.cpp

global.o: global.cpp global.h
	g++ -std=c++11 -O2 -c -o $@ global.cpp

MA.o: MA.h MPP.h MA.cpp
	g++ -std=c++11 -O2 -c -o $@ MA.cpp

clean:
	rm -f main MPP.o MA.o localsearches.o MPP_Problem.o Model.o global.o recomb.o
