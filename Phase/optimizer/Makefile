
TRANSDSPEC1DBASE=../../forwardmodel

INCLUDES = -I$(TRANSDSPEC1DBASE) $(shell gsl-config --cflags) -I$(HOME)/local/include

CXX ?= mpicxx
CXXFLAGS = -c -g -Wall -std=c++11 $(INCLUDES)

CXXFLAGS += -O3

DGGEVLIB = $(TRANSDSPEC1DBASE)/dggev/libdggev.a
SPEC1DLIB = $(TRANSDSPEC1DBASE)/spec1d/libspec1d.a

LIBS = $(SPEC1DLIB) $(DGGEVLIB) \
	$(shell gsl-config --libs) \
	-lgfortran \
	-L$(HOME)/local/lib \
	-lfftw3

ifeq ($(CXX),mpiicpc)
LIBS += -lifcore
endif

TARGETS = optimizelove \
	optimizerayleigh \
	optimizejoint

OBJS = 

all : $(TARGETS)

optimizelove: optimizelove.o $(OBJS) $(SPEC1DLIB) 
	$(CXX) -o optimizelove optimizelove.o $(OBJS) $(LIBS)

optimizerayleigh: optimizerayleigh.o $(OBJS) $(SPEC1DLIB) 
	$(CXX) -o optimizerayleigh optimizerayleigh.o $(OBJS) $(LIBS)

optimizejoint: optimizejoint.o $(OBJS) $(SPEC1DLIB) 
	$(CXX) -o optimizejoint optimizejoint.o $(OBJS) $(LIBS)

%.o : %.cpp 
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

clean :
	rm -f $(TARGETS) *.o
