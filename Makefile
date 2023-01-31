
########################################################
# modify it for different environment
ENV := mbp
#ENV := cascadia
#ENV := narval
#ENV := beluga

#compiler := intel
compiler := gnu

# modify it for different order of accuracy
pOrder := 2
########################################################

EXE  := exe_solver
EXE1 := exe_get_neigh
EXE2 := exe_part_mesh

BINDIR := ./bin
SRCDIR := ./src
OBJDIR := ./obj

EXE := $(BINDIR)/$(EXE)
EXE1 := $(BINDIR)/$(EXE1)
EXE2 := $(BINDIR)/$(EXE2)


ifeq "$(compiler)" "intel"
#FC := mpif90 -warn all -O3 -cpp -DpOrder=$(pOrder) -module $(OBJDIR)
FC := mpif90 -O3 -cpp -DpOrder=$(pOrder) -module $(OBJDIR)
endif
ifeq "$(compiler)" "gnu"
FC := mpif90 -O3 -cpp -DpOrder=$(pOrder) -J$(OBJDIR)
FC := mpif90 -Wall -O3 -cpp -DpOrder=$(pOrder) -J$(OBJDIR)
endif

#FC := $(FC) -DTW

ifeq "$(ENV)" "mbp"
######################################################
### for MacBook Pro
NETCDF_DIR := /usr/local
LDFLAGS := -lblas -llapack -L $(NETCDF_DIR)/lib -lnetcdff
#LDFLAGS := -lblas -llapack -L /usr/local/lib /usr/local/lib/libnetcdf.a /usr/local/lib/libnetcdff.a
#LDFLAGS := -L . -llapack -lrefblas 
INC := -I $(NETCDF_DIR)/include
endif

ifeq "$(ENV)" "cascadia"
######################################################
### for cascadia work station (Ubuntu 18.04)
#BLAS := /home/wzhang/spack/opt/spack/linux-ubuntu18.04-sandybridge/gcc-11.2.0/flexiblas-3.0.4-bkgkpdygjacfibzw2n4rwmdwxgc54scv
#NETCDF := /home/wzhang/spack/opt/spack/linux-ubuntu18.04-sandybridge/gcc-11.2.0/netcdf-fortran-4.5.3-4zfxvvtauqluxm5eszzs4y77s5upopw5
### for cascadia work station (Ubuntu 20.04)
BLAS := /home/wzhang/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/flexiblas-3.0.4-lsjitr2fijlzzmn73hokp4jxp325jwa2
NETCDF := /home/wzhang/spack/opt/spack/linux-ubuntu20.04-cascadelake/gcc-9.4.0/netcdf-fortran-4.5.4-cmhtaury5pjoj6devdytp3awokjpgh7s
LDFLAGS := -L $(BLAS)/lib -lflexiblas
LDFLAGS := $(LDFLAGS) -L $(NETCDF)/lib -lnetcdff
INC := -I $(NETCDF)/include
endif

ifeq "$(ENV)" "narval"
######################################################
### for compute canada (narval or beluga)
LDFLAGS := -lflexiblas -lnetcdff
INC :=
endif

ifeq "$(ENV)" "beluga"
######################################################
### for compute canada (narval or beluga)
LDFLAGS := -mkl -lnetcdff
INC :=
endif

SRC_SOLVER :=               \
       yaml_types.F90       \
       yaml.F90             \
       yaml_settings.F90    \
       mod_para.F90         \
       mod_mpi.F90          \
       mod_vtk.F90          \
       mod_math.F90         \
       mod_read.F90         \
       mod_matrix.F90       \
       mod_nodes.F90        \
       mod_types.F90        \
       mod_mesh.F90         \
       mod_jacobi.F90       \
       mod_geometry.F90     \
       mod_rotate.F90       \
       mod_eqns.F90         \
       mod_solve.F90        \
       mod_thermpress.F90   \
       mod_wave.F90         \
       mod_damp.F90         \
       mod_plastic.F90      \
       mod_source.F90       \
       mod_init_fault.F90   \
       mod_fault.F90        \
       mod_exchange.F90     \
       mod_rk.F90           \
       mod_io_fault.F90     \
       mod_io_grdsurf.F90   \
       mod_io_recv.F90      \
       mod_recv.F90         \
       mod_check.F90        \
       seis3d.F90

OBJS_SOLVER := $(foreach file,$(SRC_SOLVER),$(OBJDIR)/$(file:.F90=.o))
#SRC := $(foreach file,$(SRC),$(SRCDIR)/$(file))

.PHONY: all

all: dir solver tool

# Target
dir:
	mkdir -p bin obj
solver: $(EXE)
tool: $(EXE1) $(EXE2)

# SUFFIXES rules
.SUFFIXES:
.SUFFIXES: .F90 .o

$(EXE) : $(OBJS_SOLVER)
	$(FC) -o $@ $(OBJS_SOLVER) $(LDFLAGS)
$(EXE1) : $(OBJDIR)/get_neig.o
	$(FC) -o $@ $^ $(LDFLAGS)
$(EXE2) : $(OBJDIR)/part_mesh.o
	$(FC) -o $@ $^ $(LDFLAGS)

$(OBJDIR)/%.o : $(SRCDIR)/%.F90
	$(FC) -o $@ -c $^ $(INC)

$(OBJDIR)/get_neig.o : $(SRCDIR)/get_neigh.F90
	$(FC) -o $@ -c $^ $(INC)
$(OBJDIR)/part_mesh.o : $(SRCDIR)/part_mesh.F90
	$(FC) -o $@ -c $^ $(INC)

clean:
	rm -rf $(OBJDIR)
	mkdir $(OBJDIR)

cleanall:
	rm -rf $(BINDIR)
	rm -rf $(OBJDIR)
	mkdir $(BINDIR)
	mkdir $(OBJDIR)

