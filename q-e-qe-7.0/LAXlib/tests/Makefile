# Makefile for LAXlib

include ../../make.inc

MODFLAGS= $(MOD_FLAG)../../UtilXlib $(MOD_FLAG)../../ELPA/src  $(MOD_FLAG).. $(MOD_FLAG).

SRCS = test_diaghg_1.f90 \
       test_diaghg_2.f90 \
       test_diaghg_3.f90 \
       test_diaghg_4.f90 \
       test_diaghg_gpu_1.f90 \
       test_diaghg_gpu_2.f90 \
       test_diaghg_gpu_3.f90 \
       test_diaghg_gpu_4.f90

EXECS = $(SRCS:.f90=.x)

UTILS = tester.o mp_world.o utils.o test_io.o test_helpers.o

all: common $(EXECS)

common: $(UTILS)

%.x: %.o
	$(LD) $(LDFLAGS) $< $(UTILS) -o $@ ../libqela.a \
	../../UtilXlib/libutil.a $(QELIBS)

clean :
	- /bin/rm -f *.o *.a *.d *.i *~ *_tmp.f90 *.mod *.L *.x

# .PHONY forces execution of a rule irrespective of the presence of an
# updated file with the same name of the rule. In this way, the script 
# that generates version.f90 always runs, updating the version if you 
# execute "svn update". The update_version script takes care of not
# changing the file if the svn version did not change

.PHONY: all clean

