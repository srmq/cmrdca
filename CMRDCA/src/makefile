CXX = g++
CXXFLAGS = -O3 -march=core2 -mfpmath=sse -Wall -c -fmessage-length=0 -std=c++11 -lm -MMD -MP
OBJECTS = CMRDCA.o CMRDCAGlobal.o ConfusionMatrix.o CrispCluster.o
DEPENDS = ${OBJECTS:.o=.d}
EXEC = libCMRDCA.a

${EXEC}: ${OBJECTS}
	ar -r ${EXEC} ${OBJECTS}

clean:
	rm -f *~ *.o *.d ${EXEC}

-include ${DEPENDS}
