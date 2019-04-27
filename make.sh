rm *.o
rm *.mod

clear

module load intel

ifort -check bounds -c -O3 randgen.f params.F90 random.F90 grid.F90 ustencil.F90 lsflows.F90 vars.F90 iobelda.F90 update.F90 initconds.F90 timestep.F90
ifort -O3 main.F90 randgen.o update.o initconds.o params.o grid.o ustencil.o lsflows.o vars.o iobelda.o timestep.o random.o

module unload intel

echo 'Done'
