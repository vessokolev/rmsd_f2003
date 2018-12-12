LEVEL_1_MOD = s_count_pdb_atoms.o  s_read_pdb_atom_lines.o  s_read_pdb_coord.o  s_rmsd_center_coords.o  \
              s_rmsd_fast_calc_rmsd_and_rotation.o  s_rmsd_inner_product.o  s_rmsd_superimpose.o  \
              s_write_pdb_file.o demo.o
DEMO        = demo

ifeq "$(FC)" "f77"
FC = gfortran
endif

FCOMP = $(FC)

ifndef FCOMP
FCOMP = gfortran
endif


all: $(LEVEL_1_MOD) main

$(LEVEL_1_MOD): %.o: %.f90
	$(FCOMP) $(FFLAGS) -c $< -o $@ $(INC)

main:
	$(FCOMP) $(FFLAGS) -o $(DEMO) $(LEVEL_1_MOD)

clean:
	rm -f $(DEMO) $(LEVEL_1_MOD)

