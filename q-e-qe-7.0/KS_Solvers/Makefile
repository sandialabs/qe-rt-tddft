# Makefile for  KS_Solvers
include ../make.inc

# the following directive prevents execution of this makefile in parallel
# Problem is, libks_solvers.a depends upon all libdavid libcg etc, but if
# those dependencies are added explicitly to libks_solvers.a, the library
# is relinked every time, forcing relinking of many executables; if not, 
# parallel make may execute libks_solvers.a before all other lib* complete 
.NOTPARALLEL :

default: all

all: libdavid_rci libdavid libcg libppcg libparo librmm libdense libnewsolver libks_solvers.a ks_solver_interfaces.fh

ALLOBJS = \
Davidson_RCI/david_rci.o \
Davidson/cegterg.o \
Davidson/regterg.o \
CG/ccgdiagg.o \
CG/rcgdiagg.o \
PPCG/generic_cublas.o \
PPCG/ppcg_gamma.o \
PPCG/ppcg_k.o \
ParO/bpcg_gamma.o \
ParO/bpcg_k.o \
ParO/paro_gamma.o \
ParO/pcg_gamma.o \
ParO/paro_k.o \
ParO/pcg_k.o \
ParO/paro_gamma_new.o \
ParO/paro_k_new.o \
DENSE/rotate_HSpsi_gamma.o \
DENSE/rotate_HSpsi_k.o \
DENSE/rotate_wfc_gamma.o \
DENSE/rotate_wfc_k.o \
DENSE/rotate_driver.o \
DENSE/gram_schmidt_gamma.o \
DENSE/gram_schmidt_k.o \
DENSE/gram_schmidt_k_gpu.o \
DENSE/gram_schmidt_gamma_gpu.o \
DENSE/rotate_xpsi_gamma.o \
DENSE/rotate_xpsi_gamma_gpu.o \
DENSE/rotate_xpsi_k_gpu.o \
DENSE/rotate_xpsi_k.o \
RMM/crmmdiagg.o \
RMM/crmmdiagg_gpu.o \
RMM/rrmmdiagg.o \
RMM/rrmmdiagg_gpu.o 


# GPU-related objects
ALLOBJS += \
Davidson/cegterg_gpu.o \
Davidson/regterg_gpu.o \
DENSE/rotate_wfc_k_gpu.o \
DENSE/rotate_wfc_gamma_gpu.o \
DENSE/rotate_driver_cuf.o \
CG/rcgdiagg_gpu.o \
CG/ccgdiagg_gpu.o \
PPCG/generic_cublas.o \
PPCG/ppcg_gamma_gpu.o \
PPCG/ppcg_k_gpu.o \
ParO/bpcg_gamma_gpu.o \
ParO/bpcg_k_gpu.o \
ParO/paro_gamma_new_gpu.o \
ParO/paro_k_new_gpu.o \
DENSE/rotate_HSpsi_gamma_gpu.o \
DENSE/rotate_HSpsi_k_gpu.o

# add here other objects, e.g. ParO/*.o NewSolver/*.o

libks_solvers.a: $(ALLOBJS) 
	$(AR) $(ARFLAGS) $@ $?
	$(RANLIB) $@



libdavid_rci : 
	if test -d Davidson_RCI ; then \
	( cd Davidson_RCI ; $(MAKE) all || exit 1 ) ; fi

libdavid : 
	if test -d Davidson ; then \
        ( cd Davidson ; $(MAKE) all || exit 1 ) ; fi

libcg : 
	if test -d CG ; then \
        ( cd CG ; $(MAKE) all || exit 1 ) ; fi

libppcg : 
	if test -d PPCG ; then \
        ( cd PPCG ; $(MAKE) all || exit 1 ) ; fi

libparo : 
	if test -d ParO ; then \
        ( cd ParO ; $(MAKE) all || exit 1 ) ; fi

libdense : 
	if test -d DENSE ; then \
        ( cd DENSE ; $(MAKE) all || exit 1 ) ; fi

librmm :
	if test -d RMM ; then \
        ( cd RMM ; $(MAKE) all || exit 1 ) ; fi

libnewsolver : 
	if test -d NewSolver ; then \
        ( cd NewSolver ; $(MAKE) all || exit 1 ) ; fi

clean : 
	if test -d Davidson_RCI ; then ( cd Davidson_RCI ; $(MAKE) clean ) ; fi
	if test -d Davidson     ; then ( cd Davidson ; $(MAKE) clean ) ; fi
	if test -d CG           ; then ( cd CG ; $(MAKE) clean ) ; fi
	if test -d PPCG         ; then ( cd PPCG ; $(MAKE) clean ) ; fi
	if test -d ParO         ; then ( cd ParO ; $(MAKE) clean ) ; fi
	if test -d DENSE        ; then ( cd DENSE ; $(MAKE) clean ) ; fi
	if test -d RMM          ; then ( cd RMM ; $(MAKE) clean ) ; fi
	if test -d NewSolver    ; then ( cd NewSolver ; $(MAKE) clean ) ; fi
	rm -f *.fh 
	rm -f libks_solvers.a

 