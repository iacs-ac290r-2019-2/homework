CC := gcc
LDLIBS := -llapacke

all: main.exe

main.exe: main.o construct_matrix.o tridiag_solver.o gaussian_quadrature.o
	$(CC) -o $@ $^ $(LDLIBS)

main.o: main.c construct_matrix.h tridiag_solver.h
	$(CC) -c $< $(LDLIBS)

construct_matrix.o: construct_matrix.c construct_matrix.h
	$(CC) -c $< $(LDLIBS)

tridiag_solver.o: tridiag_solver.c tridiag_solver.h
	$(CC) -c $< $(LDLIBS)

gaussian_quadrature.o: gaussian_quadrature.c gaussian_quadrature.h
	$(CC) -c $< $(LDLIBS)

.PHONY: clean
clean:
	rm -f *.o *.exe