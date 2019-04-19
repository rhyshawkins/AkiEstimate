
BASELIB=../forwardmodel

INCLUDES = -I$(BASELIB) 

CXX ?= g++
CXXFLAGS = -c -g -Wall -std=c++11 $(INCLUDES)

CXXFLAGS += -O3

DGGEVLIB = $(BASELIB)/dggev/libdggev.a
SPEC1DLIB = $(BASELIB)/spec1d/libspec1d.a

LIBS = $(shell gsl-config --libs) $(TRANSDLIB) $(SPEC1DLIB) $(DGGEVLIB) -lgfortran

TARGETS = mkreferencerayleigh \
	mkreferencelove

OBJS = 

all : $(TARGETS)

mkreferencerayleigh: mkreferencerayleigh.o $(OBJS) $(SPEC1DLIB) $(TRANSDLIB)
	$(CXX) -o mkreferencerayleigh mkreferencerayleigh.o $(OBJS) $(LIBS)

mkreferencelove: mkreferencelove.o $(OBJS) $(SPEC1DLIB) $(TRANSDLIB)
	$(CXX) -o mkreferencelove mkreferencelove.o $(OBJS) $(LIBS)

%.o : %.cpp $(HEADERS)
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

clean :
	rm -f $(TARGETS) *.o
