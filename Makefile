
all :
	make -C forwardmodel/dggev
	make -C forwardmodel/spec1d
	make -C Reference
	make -C InitialPhase/optimizer
	make -C Phase/optimizer
	cd tutorial && pdflatex tutorial.tex 

tutorial/tutorial.pdf : tutorial/tutorial.tex
	cd tutorial && pdflatex tutorial

DISTFILES = \
	$(wildcard example_data/LoveResponse/*.txt) \
	$(wildcard example_data/RayleighResponse/*.txt) \
	$(wildcard forwardmodel/dggev/*.f) \
	forwardmodel/dggev/Makefile \
	$(wildcard forwardmodel/spec1d/*.?pp) \
	forwardmodel/spec1d/Makefile \
	$(wildcard InitialPhase/scripts/*.py) \
	$(wildcard InitialPhase/optimizer/*.?pp) \
	InitialPhase/optimizer/Makefile \
	$(wildcard Phase/optimizer/*.?pp) \
	Phase/optimizer/Makefile \
	$(wildcard Reference/*.?pp) \
	$(wildcard Reference/models/*.txt) \
	Reference/Makefile \
	$(wildcard tutorial/*.sh) \
	$(wildcard tutorial/scripts/*.py) \
	tutorial/tutorial.tex \
	tutorial/tutorial.pdf

DATE = $(shell date + "%Y%m%d%H%M")
DIR = AkiEstimate
TGZ = $(DIR).tar.gz

INSTALL = install
INSTALLFLAGS = -D

dist : tutorial/tutorial.pdf
	mkdir -p $(DIR)
	echo $(DATE) > $(DIR)/Version
	for f in Makefile $(DISTFILES); do \
	    $(INSTALL) $(INSTALLFLAGS) $$f $(DIR)/$$f ; \
	done
	tar -czf $(TGZ) $(DIR)/*
	rm -rf $(DIR)


