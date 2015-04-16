# this makefile is for LINUX machines only 

include Makefile.in

AR = ar
ITSOL2 = ./src/latools/ITSOL_2
PV = ./src/utils/pv-1.6.0
SPARSKIT2 = ./src/latools/SPARSKIT2
ITSOL2_OBJ = ${ITSOL2}/OBJ/arms2.o ${ITSOL2}/OBJ/auxill.o ${ITSOL2}/OBJ/fgmr.o \
${ITSOL2}/OBJ/indsetC.o ${ITSOL2}/OBJ/MatOps.o ${ITSOL2}/OBJ/misc.o ${ITSOL2}/OBJ/setblks.o \
${ITSOL2}/OBJ/sets.o ${ITSOL2}/OBJ/tools.o ${ITSOL2}/OBJ/arms2.o
#
DEMOSRC = example*.f90
#
FCFLAGS =  -c -g -I${ITSOL2}/INC
CCFLAGS =  -c -g -DLINUX -I${ITSOL2}/INC
#
LIB = lib
MOD = include
BIN = bin
#
DEMOSRC = demo/example_1.f90 demo/example_2.f90 \
demo/example_3.f90  demo/example_4.f90

all: src/mirana.f90 src/arms_fgmres.c
	make itsol2
	make pv
	mkdir -p ${LIB}
	mkdir -p ${MOD}
	cd ${MOD}  && \
	$(FC) ${FCFLAGS} -o ../${LIB}/mirana.o ../src/mirana.f90
	$(CC) ${CCFLAGS} -o ${LIB}/arms_fgmres.o src/arms_fgmres.c 
	make demo

demo: ${DEMOSRC} ${MOD}/mirana.mod
	mkdir -p ${BIN}
	$(FC) -I${MOD} -o ${BIN}/example_1.o demo/example_1.f90 ${LIB}/mirana.o
	$(FC) -I${MOD} -o ${BIN}/example_2.o demo/example_2.f90 ${LIB}/mirana.o
	$(FC) -I${MOD} -o ${BIN}/example_3.o demo/example_3.f90 ${LIB}/mirana.o
	$(FC) -L${ITSOL2}/LIB -o ${BIN}/example_4.o demo/example_4.f90 ${LIB}/*.o -litsol -llapack -lblas

itsol2: ${ITSOL2}
	cd ${ITSOL2} && \
	mkdir -p OBJ && \
	mkdir -p LIB && \
	make

pv: ${PV}
	cd ${PV} && \
	./configure --prefix=${MIRANA} && \
	make && \
	make install

sparskit: ${SPARSKIT2}
	cd ${SPARSKIT2} && make

doxy:
	doxygen Doxyfile

clean :
	rm -rf ${LIB} ${INCLUDE} ${BIN}
#	cd ${ITSOL2} && make cleanall
#	cd ${SPARSKIT2} && make clean

