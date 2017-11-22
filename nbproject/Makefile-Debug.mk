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
CC=gcc
CCC=g++
CXX=g++
FC=gfortran
AS=as

# Macros
CND_PLATFORM=Cygwin-Windows
CND_DLIB_EXT=dll
CND_CONF=Debug
CND_DISTDIR=dist
CND_BUILDDIR=build

# Include project Makefile
include Makefile

# Object Directory
OBJECTDIR=${CND_BUILDDIR}/${CND_CONF}/${CND_PLATFORM}

# Object Files
OBJECTFILES= \
	${OBJECTDIR}/AeroClimatology.o \
	${OBJECTDIR}/AodRetrieval.o \
	${OBJECTDIR}/AtmosphericLut.o \
	${OBJECTDIR}/InputParameter.o \
	${OBJECTDIR}/OceanReflLut.o \
	${OBJECTDIR}/S3CciL2Writer.o \
	${OBJECTDIR}/S3MetaData.o \
	${OBJECTDIR}/S3NcdfData.o \
	${OBJECTDIR}/main.o \
	${OBJECTDIR}/miscUtils.o \
	${OBJECTDIR}/tinyxml/tinystr.o \
	${OBJECTDIR}/tinyxml/tinyxml.o \
	${OBJECTDIR}/tinyxml/tinyxmlerror.o \
	${OBJECTDIR}/tinyxml/tinyxmlparser.o


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
LDLIBSOPTIONS=`pkg-config --libs netcdf-cxx4` -lm   

# Build Targets
.build-conf: ${BUILD_SUBPROJECTS}
	"${MAKE}"  -f nbproject/Makefile-${CND_CONF}.mk ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/aardvarc_slstr.exe

${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/aardvarc_slstr.exe: ${OBJECTFILES}
	${MKDIR} -p ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}
	${LINK.cc} -o ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/aardvarc_slstr ${OBJECTFILES} ${LDLIBSOPTIONS}

${OBJECTDIR}/AeroClimatology.o: AeroClimatology.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AeroClimatology.o AeroClimatology.cpp

${OBJECTDIR}/AodRetrieval.o: AodRetrieval.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AodRetrieval.o AodRetrieval.cpp

${OBJECTDIR}/AtmosphericLut.o: AtmosphericLut.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/AtmosphericLut.o AtmosphericLut.cpp

${OBJECTDIR}/InputParameter.o: InputParameter.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/InputParameter.o InputParameter.cpp

${OBJECTDIR}/OceanReflLut.o: OceanReflLut.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/OceanReflLut.o OceanReflLut.cpp

${OBJECTDIR}/S3CciL2Writer.o: S3CciL2Writer.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/S3CciL2Writer.o S3CciL2Writer.cpp

${OBJECTDIR}/S3MetaData.o: S3MetaData.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/S3MetaData.o S3MetaData.cpp

${OBJECTDIR}/S3NcdfData.o: S3NcdfData.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/S3NcdfData.o S3NcdfData.cpp

${OBJECTDIR}/main.o: main.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/main.o main.cpp

${OBJECTDIR}/miscUtils.o: miscUtils.cpp 
	${MKDIR} -p ${OBJECTDIR}
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/miscUtils.o miscUtils.cpp

${OBJECTDIR}/tinyxml/tinystr.o: tinyxml/tinystr.cpp 
	${MKDIR} -p ${OBJECTDIR}/tinyxml
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/tinyxml/tinystr.o tinyxml/tinystr.cpp

${OBJECTDIR}/tinyxml/tinyxml.o: tinyxml/tinyxml.cpp 
	${MKDIR} -p ${OBJECTDIR}/tinyxml
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/tinyxml/tinyxml.o tinyxml/tinyxml.cpp

${OBJECTDIR}/tinyxml/tinyxmlerror.o: tinyxml/tinyxmlerror.cpp 
	${MKDIR} -p ${OBJECTDIR}/tinyxml
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/tinyxml/tinyxmlerror.o tinyxml/tinyxmlerror.cpp

${OBJECTDIR}/tinyxml/tinyxmlparser.o: tinyxml/tinyxmlparser.cpp 
	${MKDIR} -p ${OBJECTDIR}/tinyxml
	${RM} "$@.d"
	$(COMPILE.cc) -g -I/cygdrive/P/C++/NR `pkg-config --cflags netcdf-cxx4`   -MMD -MP -MF "$@.d" -o ${OBJECTDIR}/tinyxml/tinyxmlparser.o tinyxml/tinyxmlparser.cpp

# Subprojects
.build-subprojects:

# Clean Targets
.clean-conf: ${CLEAN_SUBPROJECTS}
	${RM} -r ${CND_BUILDDIR}/${CND_CONF}
	${RM} ${CND_DISTDIR}/${CND_CONF}/${CND_PLATFORM}/aardvarc_slstr.exe

# Subprojects
.clean-subprojects:

# Enable dependency checking
.dep.inc: .depcheck-impl

include .dep.inc
