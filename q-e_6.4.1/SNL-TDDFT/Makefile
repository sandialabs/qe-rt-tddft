# Top level makefile for snl-tddft

all:    build

build:
	@echo "Building snl-tddft..."
	$(MAKE) -C src

clean:
	@echo "Cleaning snl-tddft..."
	if test -s src/Makefile ; then ( $(MAKE) -C src clean ); fi
	-/bin/rm -f bin/*.x

distclean:
	$(MAKE) -C src distclean
	-/bin/rm -f config.log config.status makedeps.sh
	-/bin/rm -Rf autom4te.cache
