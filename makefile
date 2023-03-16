## Makefile for COCO and Borg MOEA
CC = gcc
LDFLAGS += -lm
CCFLAGS ?= -O3 -march=native -flto -frename-registers

########################################################################
## Toplevel targets
all: experiment

clean:
	rm -f coco.o 
	rm -f experiment.o experiment 
	rm -f borg.o mt19937ar.o


########################################################################
## Programs
experiment: experiment.o coco.o borg.o mt19937ar.o
	${CC} ${CCFLAGS} -o experiment coco.o experiment.o borg.o mt19937ar.o ${LDFLAGS}  

########################################################################
## Additional dependencies
coco.o: coco.h coco.c
	${CC} -c ${CCFLAGS} -o coco.o coco.c
experiment.o: coco.h coco.c experiment.c
	${CC} -c ${CCFLAGS} -o experiment.o experiment.c
borg.o: borg.h borg.c
	${CC} -c ${CCFLAGS} -o borg.o borg.c
mt19937ar.o: mt19937ar.c mt19937ar.h
	${CC} -c ${CCFLAGS} -o mt19937ar.o mt19937ar.c
