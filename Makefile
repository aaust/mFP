MYFIT_OBJECTS=myFit.o control.o wave.o event.o likelihood.o 3j.o fitInfo.o fitInfoDict.o
MASSDEP_OBJECTS=massDep.o bw.o fitInfo.o fitInfoDict.o fitModel.o wave.o event.o
PREDICT_OBJECTS=predict.o bw.o
MOMENTS_OBJECTS=moments.o 3j.o wave.o event.o
OBJECTS=${MYFIT_OBJECTS} ${MASSDEP_OBJECTS} ${PREDICT_OBJECTS} ${MOMENTS_OBJECTS}
DEBUGFLAGS=-g -fopenmp
LDFLAGS=${DEBUGFLAGS}
CXXFLAGS:=-Wall ${DEBUGFLAGS} -O2 $(shell root-config --cflags)
LIBS:=$(shell root-config --libs --cflags) -lMinuit2 -lMathMore

all: myFit massDep predict moments
.PHONY: all

myFit: ${MYFIT_OBJECTS}
	g++ -o $@ ${LDFLAGS} ${MYFIT_OBJECTS} ${LIBS}

massDep: ${MASSDEP_OBJECTS}
	g++ -o $@ ${LDFLAGS} ${MASSDEP_OBJECTS} ${LIBS}

predict: ${PREDICT_OBJECTS}
	g++ -o $@ ${LDFLAGS} ${PREDICT_OBJECTS} ${LIBS}

moments: ${MOMENTS_OBJECTS}
	g++ -o $@ ${LDFLAGS} ${MOMENTS_OBJECTS} ${LIBS}

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
massDep.o: massDep.cc bw.h startingValue.h fitInfo.h fitModel.h
predict.o: predict.cc bw.h
moments.o: moments.cc 3j.h wave.h

control.o: control.cc control.h
wave.o: wave.cc wave.h event.h
event.o: event.cc event.h
likelihood.o: likelihood.cc likelihood.h event.h wave.h
3j.o: 3j.cc wave.h
fitInfo.o: fitInfo.cc fitInfo.h startingValue.h
fitModel.o: fitModel.cc fitModel.h bw.h
