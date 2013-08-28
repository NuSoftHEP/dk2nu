#
# $Id: GNUmakefile,v 1.1 2012-11-07 01:33:44 rhatcher Exp $
#
########################################################################

SHELL    = /bin/sh
NAME     = all
MAKEFILE = GNUmakefile

SUBDIRS = tree
### for now don't try to build GENIE interface
ifneq ($(GENIE),)
  SUBDIRS += genie
endif

all:  directories FORCE
	@echo "make all"
	for d in $(SUBDIRS) ; do (cd $$d && make all ); done
ifeq ($(GENIE),)
	@echo "GENIE not defined, skip building genie interface!!!!"
endif

clean:
	@echo "make clean"
	rm -f bin/* lib/*
	for d in $(SUBDIRS) ; do (cd $$d && make clean ); done

directories: include lib bin
	@echo "make directories"

include: 
	@echo "make include"
#   allow users to reference "dk2nu/tree/dk2nu.h" 
	ln -s .. include
#	ln -s ../dk2nu include

lib: 
	@echo "make lib"
	mkdir lib

bin: 
	@echo "make bin"
	mkdir bin

tmp:
	@echo "make tmp"
	mkdir tmp

FORCE:

########################################################################
