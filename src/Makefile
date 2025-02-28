FC:=gfortran
GFLAGS:=-O3 -cpp -std='gnu' -fopenmp
LDFLAGS:=-fopenmp
OUTPUT:=boundary_layer
SOURCES:=cal_slp_mod.f90 main.f90
OBJS:=$(patsubst %.f90, %.o, $(SOURCES))

.PHONY: all
all: $(MODS) $(OBJS)
	$(FC) $(LDFLAGS) -o $(OUTPUT) $(OBJS)

%.o: %.f90
	$(FC) $(GFLAGS) -c $< 

.PHONY: clean
clean:
	rm -f *.o *.mod $(OUTPUT) 

.PHONY: install
install:
	mv $(OUTPUT) ..
