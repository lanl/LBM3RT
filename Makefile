FC     = mpiifort
FFLAGS = -O3 -module ./mod

PROGRAM = lbm3d_mpi
SRCDIR  = src
MODDIR  = mod

SRCNAMES = \
  ptran_global.F90 \
  ptran_dbase.F90 \
  trdynmem.F90 \
  fileio.F90 \
  ptran_read.F90 \
  water_eos.F90 \
  co2eos.F90 \
  ptran_speciation.F90 \
  trgamdh.F90 \
  ptran_setbnd.F90 \
  ptran_init.F90 \
  comvar.F90 \
  setup_MPI.F90 \
  MPI_comm.F90 \
  equilf.F90 \
  init_flow.F90 \
  perind3d.F90 \
  denvel.F90 \
  streaming_flow.F90 \
  collision_flow.F90 \
  boundary_flow.F90 \
  output.F90 \
  bindex3d.F90 \
  equilg.F90 \
  solve.F90 \
  concentration.F90 \
  init_solute.F90 \
  streaming_solute.F90 \
  collision_solute.F90 \
  boundary_solute.F90 \
  expand.F90 \
  redistri_bnd.F90 \
  lbm3d.F90

SRCS = $(addprefix $(SRCDIR)/,$(SRCNAMES))

all: $(PROGRAM)

$(PROGRAM): $(MODDIR) $(SRCS)
	$(FC) $(FFLAGS) -o $(PROGRAM) $(SRCS)

$(MODDIR):
	mkdir -p $(MODDIR)
