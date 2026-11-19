SHELL =/bin/sh

DEFAULT_TARGET : general  

include ../../Makefile.def
include ../../Makefile.conf

SEARCH=

GENERALOBJ = specfunc.o rnd.o global.o 

# The immutable chemical catalogue and scalar services are application-neutral
# host C++17 code. Compile the maintained generated-tree copies through this
# existing library lane; the root archive already includes general/*.o, so no
# root Makefile or application-specific duplicate implementation is required.
# Unique object basenames prevent ar from replacing an unrelated model object.
CHEMDB_SOURCE_DIR = ../models/chemical_database
CHEMDBOBJ = ChemDbDataIO.o ChemDbRuntime.o ChemDbPhotochemistry.o \
            ChemDbElectronImpact.o ChemDbChargeExchange.o ChemDbRateServices.o
GENERALOBJ += $(CHEMDBOBJ)

ChemDb%.o: $(CHEMDB_SOURCE_DIR)/%.cpp $(wildcard $(CHEMDB_SOURCE_DIR)/*.h) input/amps_input.h
	${COMPILE.mpicxx} ${FLAGCC} $< -o $@

general : $(GENERALOBJ)
	ar -src general.a $(GENERALOBJ)
	@echo "Target 'general' is done." 
