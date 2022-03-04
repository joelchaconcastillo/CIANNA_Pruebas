main: MA.o MPP.o utils.o Model.o localsearches.o MPP_Problem.o main.cpp
	g++ -std=c++11 -O2 -o $@ main.cpp MA.o MPP.o utils.o Model.o localsearches.o MPP_Problem.o

MPP.o: MPP.h MPP.cpp
	g++ -std=c++11 -O2 -c -o $@ MPP.cpp

Model.o: MPP.h Model.cpp
	g++ -std=c++11 -O2 -c -o $@ Model.cpp

MPP_Problem.o: MPP_Problem.h MPP_Problem.cpp
	g++ -std=c++11 -O2 -c -o $@ MPP_Problem.cpp

localsearches.o: MPP.h localsearches.cpp
	g++ -std=c++11 -O2 -c -o $@ localsearches.cpp

MA.o: MA.h MPP.h MA.cpp
	g++ -std=c++11 -O2 -c -o $@ MA.cpp

utils.o: utils.cpp utils.h
	g++ -O2 -c -o $@ utils.cpp

clean:
	rm -f main MPP.o MA.o utils.o
