
FC = gfortran
FCFLAGS = -c -g -Wall -O3

CXX = g++
CXXFLAGS = -c -g -Wall

all : libdggev.a

OBJS = daxpy.o \
	dcopy.o \
	ddot.o \
	dgemm.o \
	dgemv.o \
	dgeqr2.o \
	dgeqrf.o \
	dger.o \
	dgesv.o \
	dgetrf.o \
	dgetrf2.o \
	dgetrs.o \
	dggbak.o \
	dggbal.o \
	dggev.o \
	dgghrd.o \
	dhgeqz.o \
	disnan.o \
	dlabad.o \
	dlacpy.o \
	dladiv.o \
	dlag2.o \
	dlaisnan.o \
	dlaln2.o \
	dlamch.o \
	dlange.o \
	dlanhs.o \
	dlapy2.o \
	dlapy3.o \
	dlarf.o \
	dlarfb.o \
	dlarfg.o \
	dlarft.o \
	dlartg.o \
	dlascl.o \
	dlaset.o \
	dlassq.o \
	dlasv2.o \
	dlaswp.o \
	dlasyf.o \
	dnrm2.o \
	dorg2r.o \
	dorgqr.o \
	dorm2r.o \
	dormqr.o \
	drot.o \
	dscal.o \
	dswap.o \
	dtgevc.o \
	dtrmm.o \
	dtrmv.o \
	dtrsm.o \
	dsyconv.o \
	dsyr.o \
	dsysv.o \
	dsytf2.o \
	dsytrf.o \
	dsytrs.o \
	dsytrs2.o \
	idamax.o \
	ieeeck.o \
	iladlc.o \
	iladlr.o \
	ilaenv.o \
	iparam2stage.o \
	iparmq.o \
	lsame.o \
	xerbla.o

SOURCES = daxpy.f \
	dcopy.f \
	ddot.f \
	dgemm.f \
	dgemv.f \
	dgeqr2.f \
	dgeqrf.f \
	dger.f \
	dgesv.f \
	dgetrf.f \
	dgetrf2.f \
	dgetrs.f \
	dggbak.f \
	dggbal.f \
	dggev.f \
	dgghrd.f \
	dhgeqz.f \
	disnan.f \
	dlabad.f \
	dlacpy.f \
	dladiv.f \
	dlag2.f \
	dlaisnan.f \
	dlaln2.f \
	dlamch.f \
	dlange.f \
	dlanhs.f \
	dlapy2.f \
	dlapy3.f \
	dlarf.f \
	dlarfb.f \
	dlarfg.f \
	dlarft.f \
	dlartg.f \
	dlascl.f \
	dlaset.f \
	dlassq.f \
	dlasv2.f \
	dlaswp.f \
	dlasyf.f \
	dnrm2.f \
	dorg2r.f \
	dorgqr.f \
	dorm2r.f \
	dormqr.f \
	drot.f \
	dscal.f \
	dswap.f \
	dtgevc.f \
	dtrmm.f \
	dtrmv.f \
	dtrsm.f \
	dsyconv.f \
	dsyr.f \
	dsysv.f \
	dsytf2.f \
	dsytrf.f \
	dsytrs.f \
	dsytrs2.f \
	idamax.f \
	ieeeck.f \
	iladlc.f \
	iladlr.f \
	ilaenv.f \
	iparam2stage.f \
	iparmq.f \
	lsame.f \
	xerbla.f

libdggev.a : $(OBJS)
	ar -r $@ $(OBJS)

%.o : %.f
	$(FC) $(FCFLAGS) -o $*.o $*.f


linktest : libdggev.a linktest.o
	$(CXX) -o $@ linktest.o libdggev.a -lgfortran

linktest.o : linktest.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< 

INSTALL = install
INSTALLFLAGS = -D

DATE = $(shell date +"%Y%m%d%H%M")
DIR = dggev
TGZ = $(DIR).tar.gz

dist :
	mkdir -p $(DIR)
	echo $(DATA) > $(DIR)/Version
	for f in Makefile $(SOURCES); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)

clean :
	rm -f libdggev.a *.o
