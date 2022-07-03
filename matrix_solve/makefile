CC=mpicc

CFLAGS=-O0 -fPIC -mieee -mftz
INCLUDE=-I.

LIB=libmatrix_solve.a
EXE=matrix_solve

all: $(LIB) $(EXE)

$(EXE): main.o host.o slave.o
	mpicxx -mhybrid -o $(EXE) $^ -L. -lmatrix_solve

main.o:	main.cpp
	$(CC) -mhost $(CFLAGS) -fpermissive $(INCLUDE) -c $< -o $@ 

host.o: host.cpp
	$(CC) -mhost $(CFLAGS) -fpermissive $(INCLUDE) -c $< -o $@	

slave.o: slave.c
	$(CC) -mslave -msimd $(CFLAGS) -c $< -o $@

clean:
	$(RM) $(EXE) *.o
