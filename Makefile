OBJECTS=myFit.o control.o wave.o event.o likelihood.o 3j.o fitInfo.o fitInfoDict.o
DEBUGFLAGS=-g -fopenmp
LDFLAGS=${DEBUGFLAGS}
CXXFLAGS:=-Wall ${DEBUGFLAGS} -O2 $(shell root-config --cflags)
LIBS:=$(shell root-config --libs --cflags) -lMinuit2 -lMathMore

all: myFit
.PHONY: all

myFit: ${OBJECTS}
	g++ -o $@ ${LDFLAGS} ${OBJECTS} ${LIBS}

.PHONY: clean
clean:
	rm -f ${OBJECTS}
	rm -f fitInfoDict.{cc,h}
	rm -f myFit

.cc.o:
	g++ -c ${CXXFLAGS} $< -o $@

fitInfoDict.o: fitInfoDict.cc
fitInfoDict.cc: fitInfo.h startingValue.h
	rootcint -f $@ -c fitInfo.h+

myFit.o: myFit.cc control.h wave.h likelihood.h gHist.h startingValue.h fitInfo.h
control.o: control.cc control.h
wave.o: wave.cc wave.h event.h
event.o: event.cc event.h
likelihood.o: likelihood.cc likelihood.h event.h wave.h
3j.o: 3j.cc wave.h
fitInfo.o: fitInfo.cc fitInfo.h startingValue.h
