CC := gcc

all: main.exe

main.exe: main.o tridiag_solver.o
	$(CC) -o $@ $^ -llapacke

main.o: main.c tridiag_solver.h
	$(CC) -c $< -llapacke

tridiag_solver.o: tridiag_solver.c tridiag_solver.h
	$(CC) -c $< -llapacke

.PHONY: clean
clean:
	rm -f *.o *.exe