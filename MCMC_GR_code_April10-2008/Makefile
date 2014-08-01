
CC=gcc
CFLAGS=-c -Wall

all:hash

hash: HASHcode.o mcmcsampling.o fragmatrix.o
	$(CC)  HASHcode.o fragmatrix.o mcmcsampling.o -o HASH -lm 

fragmatrix.o: fragmatrix.h fragmatrix.c common.h
	$(CC) -c fragmatrix.c

mcmcsampling.o: mcmcsampling.h mcmcsampling.c common.h
	$(CC) -c mcmcsampling.c

HASHcode.o: HASHcode.c mcmcsampling.h mcmcsampling.c common.h fragmatrix.h fragmatrix.c
	$(CC) -c HASHcode.c

clean: 
	rm -f fragmatrix.o mcmcsampling.o HASHcode.o HASH 




