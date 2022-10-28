gcc_options = -std=c++17 -Wall --pedantic-errors -DMKL_ILP64  -I"${MKLROOT}/include" -g
l_b = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

program : DNS_HamiltonianTest.o DNS_Hamiltonian.o Jset.o EIGEN.o
	g++ -o $@ $^ $(l_b)

DNS_HamiltonianTest.o : DNS_HamiltonianTest.cpp
	g++ -c $< $(l_b)

DNS_Hamiltonian.o : DNS_Hamiltonian/DNS_Hamiltonian.cpp
	g++ -c $< $(l_b)

Jset.o : Jset/Jset.cpp
	g++ -c $< $(l_b)
	
EIGEN.o : EIGEN/EIGEN.cpp
	g++ -c $< $(l_b)

run : program
	./program

clean:
	rm -f ./program

.PHONY : run clean
