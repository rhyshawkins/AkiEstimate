
CXX ?= g++
CXXFLAGS = -c -g -Wall -std=c++11 $(INCLUDES)

TARGETS = libspec1d.a

LIBS = $(shell gsl-config --libs)

AR = ar
ARFLAGS = -r

INSTALL = install
INSTALLFLAGS = -D

OBJS = ak135.o \
	iasp91.o \
	logging.o \
	simplified_ak135.o

all : $(TARGETS)

libspec1d.a : $(OBJS)
	$(AR) $(ARFLAGS) $@ $(OBJS)	

%.o : %.cpp
	$(CXX) $(CXXFLAGS) -o $*.o $*.cpp

SRCS = Makefile \
	ak135.hpp \
	cell.hpp \
	density.hpp \
	eigenroots.hpp \
	empiricalmodel.hpp \
	encodedecode.hpp \
	fixedboundary.hpp \
	generalisedeigenproblem.hpp \
	iasp91.hpp \
	isotropicrhovpvs.hpp \
	isotropicrhovpvshalfspace.hpp \
	isotropicrhovs.hpp \
	isotropicrhovshalfspace.hpp \
	isotropicvs.hpp \
	isotropicvshalfspace.hpp \
	laguerre.hpp \
	laguerrequadrature.hpp \
	legendre.hpp \
	lobatto.hpp \
	lobattoprojection.hpp \
	lobattoquadrature.hpp \
	logging.hpp \
	loveanalytical.hpp \
	lovematrices.hpp \
	mesh.hpp \
	model.hpp \
	modelhistory.hpp \
	modelinterface.hpp \
	modelloader.hpp \
	parameterset.hpp \
	polynomial.hpp \
	rayleighmatrices.hpp \
	regression.hpp \
	simplified_ak135.hpp \
	spec1dmatrix.hpp \
	ak135.cpp \
	iasp91.cpp \
	logging.cpp \
	simplified_ak135.cpp



DATE = $(shell date +"%Y%m%d%H%M")
DIR = spec1d
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in $(SRCS) $(EXTRADIST); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean :
	rm -f $(TARGETS) *.o

