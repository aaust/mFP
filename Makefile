OBJECTS=myFit.o control.o wave.o event.o likelihood.o 3j.o
LDFLAGS=-g
CXXFLAGS=-Wall -g -O2 `root-config --cflags`
LIBS=`root-config --libs --cflags` -lMinuit2 -lMathMore

myFit: ${OBJECTS}
	g++ -o $@ ${LDFLAGS} ${OBJECTS} ${LIBS}

clean:
	rm -f ${OBJECTS}
	rm -f myFit

.cc.o:
	g++ -c ${CXXFLAGS} $< -o $@

myFit.o: myFit.cc control.h wave.h
control.o: control.cc control.h
wave.o: wave.cc wave.h event.h
event.o: event.cc event.h
likelihood.o: likelihood.cc likelihood.h event.h wave.h
3j.0: 3j.cc wave.h
