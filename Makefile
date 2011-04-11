myFit: myFit.cc
	g++ -Wall -o myFit -O2 myFit.cc `root-config --libs --cflags` -lMinuit2 -lMathMore
