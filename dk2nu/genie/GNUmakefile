#
# $Id: GNUmakefile,v 1.2 2012-11-21 04:48:40 rhatcher Exp $
#
########################################################################

SHELL    = /bin/sh
NAME     = all
MAKEFILE = GNUmakefile

PACKAGE  = dk2nuGenie
LIBDEPS  = -ldk2nuTree -lPhysics   # TVector3, TLorentzVector

all: lib
	@echo "make all $(PACKAGE)"

include $(DK2NU)/make.include

lib: $(SHAREDLIB) 
	@echo "make lib requires $(SHAREDLIB)"

########################################################################
