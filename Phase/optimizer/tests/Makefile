
TRANSDSPEC1DBASE=../../../forwardmodel

INCLUDES = -I../ -I$(TRANSDSPEC1DBASE) $(shell gsl-config --cflags) -I$(HOME)/local/include

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

TARGETS = test_joint_skip

OBJS = 

all : $(TARGETS)

test_joint_skip: test_joint_skip.o $(OBJS) $(SPEC1DLIB) 
	$(CXX) -o test_joint_skip test_joint_skip.o $(OBJS) $(LIBS)

%.o : %.cpp 
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

clean :
	rm -f $(TARGETS) *.o
