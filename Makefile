CC=mpicc
CFLAGS=-O2

tspgen:
	$(CC) -O3 -o tspgen *.c -lm

.PHONY: clean
clean:
	rm -f tspgen

.PHONY: all
all:
	tspgen
