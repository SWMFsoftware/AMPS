SHELL=/bin/sh

DEFAULT_TARGET : amps

# These definitions may be overwritten by Makefile.def
SOURCES=src
WSD=srcTemp
SPICE=nospice
OPENMP=off

#extra compiler options can be defined in Makefile.local
EXTRACOMPILEROPTIONS=

#extra linker options specific for fortran and c++ linker
EXTRALINKEROPTIONS_F=
EXTRALINKEROPTIONS_CPP=

#Compiling with the CCMC's Kameleon
KAMELEON=nokameleon
BOOST=noboost
BATL=nobatl
TESTMODE=off
INTERFACE=off

#Link the SWMF' shared library
LINK_SWMF_SHARED_LIB=off

include Makefile.conf
include Makefile.def


#include the local $(MAKE)file (defined the AMPS' compiling variables)  
include Makefile.local

#the default value of the c++ compiler flags
SEARCH_C=-DMPI_ON -LANG:std -I${CWD}/${WSD}/pic -I${CWD}/${WSD}/main  -I${CWD}/${WSD}/meshAMR -I${CWD}/${WSD}/interface -I${CWD}/${WSD}/general -I${CWD}/${WSD}/models/electron_impact -I${CWD}/${WSD}/models/sputtering -I${CWD}/${WSD}/models/dust -I${CWD}/${WSD}/models/charge_exchange -I${CWD}/${WSD}/models/photolytic_reactions -I${CWD}/${WSD}/species -I${CWD}/${WSD}/models/exosphere -I${CWD}/${WSD}/models/surface -I${SPICE}/include -I${BOOST}/include -I${KAMELEON}/src -I${CWD}/utility/PostProcess -I${CWD}/share/Library/src  -I${CWD}

SEARCH_C+=${EXTRACOMPILEROPTIONS}

SEARCH_C_GENERAL=

#define the "compile kameleon' flag only when KAMELEON is used (to exclude including of the KAMELEON headers on machimes where KAMELEON is not installed) 
ifneq ($(KAMELEON),nokameleon)   
	SEARCH_C+=-D _PIC_COMPILE__KAMELEON_ 
else
	SEARCH_C+=-D _NO_KAMELEON_CALLS_
endif

#define _NO_SPICE_CALLS_ when SPICE library is not intended to be linked
ifeq ($(SPICE),nospice)
	SEARCH_C+=-D _NO_SPICE_CALLS_
endif


#the additional argument string for the fortran compiler
SEARCH_F=
#-fdefault-real-8 

# These definitions are inherited from Makefile.def and Makefile.conf
CC=${COMPILE.mpicxx}
CWD=${PTDIR}
AMPSLINKER=${CC}

AMPSLINKLIB= 

EXTRALINKEROPTIONS=

ifeq ($(LINK_SWMF_SHARED_LIB),on)
	AMPSLINKER=${LINK.f90}
	AMPSLINKLIB+=./share/lib/libSHARE.a
endif

#include BATL-related libraries for linking
ifneq ($(BATL),nobatl)	
	AMPSLINKLIB+=${BATL}/lib/libREADAMR.a
	AMPSLINKLIB+=${BATL}/lib/libSHARE.a
	AMPSLINKER=${LINK.f90}

ifeq ($(COMPILE.f90),gfortran)  
	SEARCH_F+= -J${BATL}/share/include 
else ifeq ($(COMPILE.f90),mpif90)
	SEARCH_F+= -I${BATL}/share/include
else ifeq ($(COMPILE.f90),pgf90)
	SEARCH_F+= -module ${BATL}/share/include
else ifeq  ($(COMPILE.f90),ifort)
	SEARCH_F+= -module ${BATL}/share/include
else ifeq  ($(COMPILE.f90),nagfor)
	SEARCH_F+= -I${BATL}/share/include
endif

endif

# include interface with external FORTRAN subroutines
ifeq ($(INTERFACE),on)
	AMPSLINKER=${LINK.f90}
	AMPSLINKLIB+=${WSD}/interface/interface.a	
endif 

#when OpenMP is used add the appropriate compiler flag and library
ifeq ($(OPENMP),on)
	SEARCH_C+=-fopenmp
	SEARCH_C_GENERAL+=-fopenmp
	AMPSLINKLIB+=-fopenmp
endif


# when linking mixed C/C++ and FORTRAN code mpif90 is used for linking
# certain compilers (Intel, PGI) require an additional flag
# if main() subroutine is written in C/C++

ifeq (${AMPSLINKER},${LINK.f90})

ifeq ($(COMPILE.f90),${LINK.f90})
	AMPSLINKER=mpif90  
endif

ifeq ($(COMPILE.f90),pgf90)
	AMPSLINKLIB+= -Mnomain 
else ifeq  ($(COMPILE.f90),ifort)
	AMPSLINKLIB+= -nofor-main  
else ifeq  ($(COMPILE.f90),gfortran)
	AMPSLINKLIB+= 
endif
endif

#add liker options that are apecific cor c++ cortran linkers 
ifeq (${AMPSLINKER},mpif90) 
  AMPSLINKLIB+=${EXTRALINKEROPTIONS_F} 
else 
  AMPSLINKLIB+=${EXTRALINKEROPTIONS_CPP} 
endif 


#add KAMELEON to the list of the linked libraries
ifneq ($(KAMELEON),nokameleon)
	AMPSLINKLIB+=/Users/ccmc/miniconda2/envs/kameleon/lib/ccmc/libccmc.dylib  
endif

#add SPICE to the list of the linked libraries 
ifneq ($(SPICE),nospice)
	AMPSLINKLIB+=${SPICE}/lib/cspice.a
endif

#set the nightly test flag
ifeq ($(TESTMODE),on)  
	SEARCH_C+=-D _PIC_NIGHTLY_TEST_MODE_=_PIC_MODE_ON_
endif

install:
	rm -f output
	ln -s data/output output
	@echo "AMPS installed"


distclean:
	./Config.pl -uninstall

allclean: clean
	rm -rf main srcTemp *.input* amps Makefile.local Makefile.test \
		.amps.conf .general.conf
	rm -f output

rundir:
	mkdir -p ${RUNDIR}/PT
	cd ${RUNDIR}/PT; mkdir restartIN restartOUT plots


EXE=amps

LIB_AMPS = ${WSD}/libAMPS.a

clean:
	rm -rf ${LIB_AMPS} ${WSD}
	@(if [ -d srcInterface ]; then cd srcInterface; $(MAKE) clean; fi);

tar:
	cd ../pic-tower/sources/general; rm -f *.o *.a
	cd ../pic-tower/sources/dsmc; rm -f *.o *.a
	tar -cvf sources.tar sources

${WSD}:
	./ampsConfig.pl -input ${InputFileAMPS} -no-compile 
	./utility/CheckMacro.pl ${WSD} -in-place

LIB: 
	@(if [ -d ${WSD} ]; then rm -rf ${WSD}; fi)
	$(MAKE) ${WSD}
	$(MAKE) LIB_after_build
	(if [ "$(STANDALONE)" == "NO" ]; then \
		cd srcInterface; $(MAKE) LIB SEARCH_C="${SEARCH_C}"; fi)

LIB_after_build: 
ifeq ($(INTERFACE),on)
	cd ${WSD}/interface; $(MAKE) SEARCH_C="${SEARCH_C}" SEARCH="${SEARCH_F}" 
endif
	cd ${WSD}/general;                     $(MAKE) SEARCH_C="${SEARCH_C_GENERAL}" 
	cd ${WSD}/meshAMR;                     $(MAKE) SEARCH_C="${SEARCH_C}" 
	cd ${WSD}/pic;                         $(MAKE) SEARCH_C="${SEARCH_C}" SEARCH="${SEARCH_F}" 
	cd ${WSD}/species;                     $(MAKE) SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/exosphere;            $(MAKE) SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/surface;              $(MAKE) SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/electron_impact;      $(MAKE) SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/sputtering;           $(MAKE) SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/dust;                 $(MAKE) SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/charge_exchange;      $(MAKE) SEARCH_C="${SEARCH_C}"
	cd ${WSD}/models/photolytic_reactions; $(MAKE) SEARCH_C="${SEARCH_C}" 
	cd ${WSD}/main; $(MAKE) SEARCH_C="${SEARCH_C}"
	cp -f ${WSD}/main/mainlib.a ${WSD}/libAMPS.a
ifeq ($(SPICE),nospice)
	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o models/exosphere/*.o models/surface/*.o models/electron_impact/*.o models/sputtering/*.o models/dust/*.o models/charge_exchange/*.o models/photolytic_reactions/*.o
else
	rm -rf ${WSD}/tmpSPICE
	mkdir ${WSD}/tmpSPICE
	cp ${SPICE}/lib/cspice.a ${WSD}/tmpSPICE
	cd ${WSD}/tmpSPICE; ar -x cspice.a
	cd ${WSD}; ${AR} libAMPS.a general/*.o meshAMR/*.o pic/*.o species/*.o models/exosphere/*.o models/surface/*.o models/electron_impact/*.o models/sputtering/*.o models/dust/*.o models/charge_exchange/*.o models/photolytic_reactions/*.o tmpSPICE/*.o
endif

.PHONY: amps
amps:
	@(if [ -d ${WSD} ]; then rm -rf ${WSD}; fi)
	$(MAKE) $(WSD)
	$(MAKE) amps_after_build

amps_after_build: LIB_after_build 
	@rm -f amps
	cd ${WSD}/main; $(MAKE) amps SEARCH_C="${SEARCH_C}"
	$(MAKE) amps_link

amps_link:
	${AMPSLINKER} -o amps srcTemp/main/main.a srcTemp/libAMPS.a \
		${CPPLIB} ${AMPSLINKLIB} ${EXTRALINKEROPTIONS}

.PHONY: test
test:
	echo "Use make test_all to run the AMPS tests"
