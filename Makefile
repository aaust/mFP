OBJECTS=myFit.o control.o wave.o event.o likelihood.o
LDFLAGS=-g
CXXFLAGS=-Wall -g -O2
LIBS=`root-config --libs --cflags` -lMinuit2 -lMathMore

myFit: ${OBJECTS}
	g++ -o $@ ${LDFLAGS} ${OBJECTS} ${LIBS}

.cc.o:
	g++ -c ${CXXFLAGS} $< -o $@

myFit.o: myFit.cc control.h wave.h
control.o: control.cc control.h
wave.o: wave.cc wave.h event.h
event.o: event.cc event.h
likelihood.o: likelihood.cc likelihood.h event.h wave.h
