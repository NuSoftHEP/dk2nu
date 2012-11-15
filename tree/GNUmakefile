#
# $Id: GNUmakefile,v 1.3 2012-11-15 09:09:27 rhatcher Exp $
#
########################################################################

SHELL    = /bin/sh
NAME     = all
MAKEFILE = GNUmakefile

PACKAGE  = dk2nuTree
LIBDEPS  = -lPhysics    # TVector3

all:  FORCE lib
	@echo "make all $(PACKAGE)"

include $(DK2NU)/make.include

lib: $(SHAREDLIB) 
	@echo "make lib requires $(SHAREDLIB)"

FORCE:
	@echo "FORCE SHAREDLIB:  $(SHAREDLIB)"
	@echo "FORCE dictionary: $(DICTIONARY)"
	@echo "FORCE objects:    $(objects)"
	@echo "FORCE headers:    $(headers)"
	@echo "FORCE sources:    $(sources)"
	@echo "FORCE depfiles:   $(depfiles)"
	@echo "FORCE INCLUDE:    $(INCLUDE)"
	@echo "FORCE RCOREARCH:  $(RCOREARCH)"
	@echo "FORCE IS32:       \"$(IS32)\""
	@echo "FORCE CXXFLAGS:   $(CXXFLAGS)"
	@echo "FORCE LDFLAGS:    $(LDFLAGS)"
########################################################################
