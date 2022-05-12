all: ga

ga: curve.cpp curve.h main.cpp
	g++ -o bezier curve.cpp main.cpp -lm -lGL -lGLU -lglut -lGLEW
	
run: ga
	./bezier < rr.in > rr.out

clean:
	rm ga
