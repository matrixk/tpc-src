CC  = gcc
CXX = g++
INCLUDE =
CFLAGS = -Wall -O2
LIBS = -lm
EXE = geant2root 2

all: $(EXE)

# INCLUDE += $(shell gsl-config --cflags)
# LIBS    += $(shell gsl-config --libs)

ROOTCFLAGS      = $(shell root-config --cflags)
ROOTLIBS        = $(shell root-config --libs)
ROOTGLIBS       = $(shell root-config --glibs)

CFLAGS += $(ROOTCFLAGS)
LIBS   += $(ROOTLIBS)
LIBS   += $(ROOTGLIBS)

#DrawSingleGamma:DrawSingleGamma.C
#DrawElectronTrack:DrawElectronTrack.C
# $(EXE): $(EXE).C ../lib/lib.h


#DrawGaussElectronPositin: DrawGaussElectronPosition.C
#DrawGamma:DrawGamma.C
#DrawGaussLine:DrawGaussLine.C
#DrawGammaMultiScatter:DrawGammaMultiScatter.C
geant2root: geant2root.C ../lib/lib.h
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBS) $(LDFLAGS) -o $@

2: 2.C ../lib/lib.h
	$(CXX) $(CFLAGS) $(INCLUDE) $< $(LIBS) $(LDFLAGS) -o $@

LinkDef.h: Makefile
	echo '#pragma link C++ class vector<int>+;' > $@
	echo '#pragma link C++ class vector<vector<float> >+;' >> $@
Vector.h: Makefile
	echo '#include <vector>' > $@
VectorDict.cpp: Vector.h LinkDef.h
	rootcint -f $@ -c $^
VectorDict.o: VectorDict.cpp VectorDict.h
	$(CC) $(CFLAGS) $(INCLUDE) $< -c -o $@

clean:
	rm -f *.o $(EXE)
