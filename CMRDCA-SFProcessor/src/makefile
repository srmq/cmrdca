USER_OBJDIR = ../../CMRDCA/src
CXX = g++
CXXFLAGS = -I"${USER_OBJDIR}" -O3 -march=core2 -mfpmath=sse -fopenmp -Wall -c -fmessage-length=0 -std=c++11 -lm -MMD -MP
LDFLAGS = -lm -lgomp
OBJECTS = CMRDCA-SFProcessor.o
USER_OBJS = ${USER_OBJDIR}/libCMRDCA.a
DEPENDS = ${OBJECTS:.o=.d}
EXEC = CMRDCA-SFProcessor

${EXEC}: ${OBJECTS} ${USER_OBJS}
	${CXX} ${OBJECTS} ${USER_OBJS} ${LDFLAGS} -o ${EXEC}

${USER_OBJS}:
	$(MAKE) -C ${USER_OBJDIR}

clean:
	rm -f *~ *.o *.d ${EXEC} ${USER_OBJDIR}/*.o ${USER_OBJDIR}/*.d ${USER_OBJS} ${USER_OBJDIR}/*~

-include ${DEPENDS}
