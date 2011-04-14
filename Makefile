OBJECTS=myFit.o control.o
LDFLAGS=-g
CXXFLAGS=-Wall -g -O2
LIBS=`root-config --libs --cflags` -lMinuit2 -lMathMore

myFit: ${OBJECTS}
	g++ -o $@ ${LDFLAGS} ${OBJECTS} ${LIBS}

.cc.o:
	g++ -c ${CXXFLAGS} $< -o $@

myFit.o: myFit.cc control.h
control.o: control.cc control.h
