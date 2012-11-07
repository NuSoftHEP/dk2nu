#
# $Id: GNUmakefile,v 1.1 2012-11-07 01:36:55 rhatcher Exp $
#
########################################################################

SHELL    = /bin/sh
NAME     = all
MAKEFILE = GNUmakefile

PACKAGE  = dk2nuGenie
LIBDEPS  = -ldk2nuTree -lEG   # TVector3

all: lib
	@echo "make all $(PACKAGE)"

include $(DK2NU)/make.include

lib: $(SHAREDLIB) 
	@echo "make lib requires $(SHAREDLIB)"

########################################################################
