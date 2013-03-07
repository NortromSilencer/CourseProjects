#
# Generated Makefile - do not edit!
#
# Edit the Makefile in the project folder instead (../Makefile). Each target
# has a -pre and a -post target defined where you can add customized code.
#
# This makefile implements configuration specific macros and targets.


# Environment
MKDIR=mkdir
CP=cp
GREP=grep
NM=nm
CCADMIN=CCadmin
RANLIB=ranlib
CC=gcc.exe
CCC=g++.exe -ffree-form
CXX=g++.exe -ffree-form
FC=gfortran.exe
AS=as.exe

# Macros
CND_PLATFORM=Cygwin-Windows
CND_DLIB_EXT=dll
CND_CONF=Release
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/sigcan.o \
	${OBJECTDIR}/etkernel.o \
	${OBJECTDIR}/efi.o \
	${OBJECTDIR}/ggubfs.o \
	${OBJECTDIR}/lambda.o \
	${OBJECTDIR}/sigsample.o \
	${OBJECTDIR}/energy.o \
	${OBJECTDIR}/angles.o \
	${OBJECTDIR}/valinterp.o \
	${OBJECTDIR}/locate.o \
	${OBJECTDIR}/rusroulette.o \
	${OBJECTDIR}/angdist.o \
	${OBJECTDIR}/intsqw.o \
	${OBJECTDIR}/newpos.o \
	${OBJECTDIR}/foo.o


# C Compiler Flags
CFLAGS=

# CC Compiler Flags
CCFLAGS=
CXXFLAGS=

# Fortran Compiler Flags
FFLAGS=

# Assembler Flags
ASFLAGS=

# Link Libraries and Options
LDLIBSOPTIONS=

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mcmarilps.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mcmarilps.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.f} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mcmarilps ${OBJECTFILES} ${LDLIBSOPTIONS} 

${OBJECTDIR}/sigcan.o: sigcan.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/sigcan.o sigcan.f90

${OBJECTDIR}/etkernel.o: etkernel.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/etkernel.o etkernel.f90

${OBJECTDIR}/efi.o: efi.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/efi.o efi.f90

${OBJECTDIR}/ggubfs.o: ggubfs.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/ggubfs.o ggubfs.f90

${OBJECTDIR}/lambda.o: lambda.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/lambda.o lambda.f90

${OBJECTDIR}/sigsample.o: sigsample.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/sigsample.o sigsample.f90

${OBJECTDIR}/energy.o: energy.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/energy.o energy.f90

${OBJECTDIR}/angles.o: angles.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/angles.o angles.f90

${OBJECTDIR}/valinterp.o: valinterp.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/valinterp.o valinterp.f90

${OBJECTDIR}/locate.o: locate.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/locate.o locate.f90

${OBJECTDIR}/rusroulette.o: rusroulette.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/rusroulette.o rusroulette.f90

${OBJECTDIR}/angdist.o: angdist.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/angdist.o angdist.f90

${OBJECTDIR}/intsqw.o: intsqw.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/intsqw.o intsqw.f90

${OBJECTDIR}/newpos.o: newpos.f90 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/newpos.o newpos.f90

${OBJECTDIR}/foo.o: foo.f 
	${MKDIR} -p ${OBJECTDIR}
	$(COMPILE.f) -O2 -o ${OBJECTDIR}/foo.o foo.f

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/mcmarilps.exe
	${RM} *.mod

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
