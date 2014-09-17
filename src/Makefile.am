EXTRA_DIST = example
# DIST_SUBDIRS = $(SUBDIRS) example

DEFS = $(REVISION)

if WANT_COLTREE
bin_PROGRAMS = coltree
PROG = COLTREE
else
if WANT_EVOLVE
bin_PROGRAMS = evolve
PROG = EVOLVE
else
bin_PROGRAMS = fitmodel
PROG = FITMODEL
endif
endif

if WANT_EVOLVE
evolve_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
models.c  models.h\
free.c  free.h\
options.c  options.h\
eigen.c  eigen.h\
eigenmb.c  eigenmb.h\
draw.c  draw.h
evolve_LDADD = -lm 
else
if WANT_COLTREE
coltree_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
models.c  models.h\
free.c  free.h\
options.c  options.h\
eigen.c  eigen.h\
eigenmb.c  eigenmb.h\
draw.c  draw.h
coltree_LDADD = -lm 
else
fitmodel_SOURCES = main.c \
utilities.c  utilities.h\
optimiz.c  optimiz.h\
lk.c  lk.h\
models.c  models.h\
free.c  free.h\
options.c  options.h\
eigen.c  eigen.h\
eigenmb.c  eigenmb.h\
draw.c  draw.h
fitmodel_LDADD = -lm 
endif
endif


all-am:	intro $(bin_PROGRAMS)
	@echo ""
	@echo "Done."

intro:	
	@echo ""
	@echo ""
	@echo ".: Building [$(bin_PROGRAMS)]. Version $(VERSION) :."
	@echo ""
	@echo ""

