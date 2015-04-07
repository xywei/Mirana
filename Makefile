# this makefile is for LINUX machines only 

AR = ar
ITSOL2 = ./src/latools/ITSOL_2
SPARSKIT2 = ./src/latools/SPARSKIT2
#
DEMOSRC = example*.f90
#
FC      =  ifort
FCFLAGS =  -c -g -I${ITSOL2}/INC
CC      =  icc
CCFLAGS =  -c -g -DLINUX -I${ITSOL2}/INC
#
LIB = lib
MOD = include
BIN = bin
#
DEMOSRC = demo/example_1.f90 demo/example_2.f90 \
demo/example_3.f90  demo/example_4.f90

all: src/mirana.f90 src/arms_fgmres.c ${DEMOSRC} ${ITSOL2}/LIB/libitsol.a
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
	$(FC) -I${MOD} -I${ITSOL2}/INC -L${ITSOL2}/OBJ -o ${BIN}/example_4.o demo/example_4.f90 ${LIB}/*.o ${ITSOL2}/LIB/libitsol.a -llapack -lblas

clean :
	rm -rf ${LIB} ${INCLUDE} ${BIN}

