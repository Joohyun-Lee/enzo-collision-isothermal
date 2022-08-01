#include <stdio.h>
void auto_show_flags(FILE *fp) {
   fprintf (fp,"\n");
   fprintf (fp,"CPP = /usr/bin/cpp\n");
   fprintf (fp,"CC  = /home/jhl1862/library/openmpi-4/bin/mpicc\n");
   fprintf (fp,"CXX = /home/jhl1862/library/openmpi-4/bin/mpic++\n");
   fprintf (fp,"FC  = /usr/bin/gfortran\n");
   fprintf (fp,"F90 = /usr/bin/gfortran\n");
   fprintf (fp,"LD  = /home/jhl1862/library/openmpi-4/bin/mpic++\n");
   fprintf (fp,"\n");
   fprintf (fp,"DEFINES = -DLINUX -DH5_USE_16_API   -D__max_subgrids=100000 -D__max_baryons=30 -D__max_cpu_per_node=48 -D__memory_pool_size=100000 -DINITS64 -DLARGE_INTS -DCONFIG_PINT_8 -DIO_32   -DNEW_PROBLEM_TYPES -DUSE_MPI   -DCONFIG_PFLOAT_8 -DCONFIG_BFLOAT_8  -DUSE_HDF5_GROUPS   -DTRANSFER   -DNEW_GRID_IO -DFAST_SIB      -DENZO_PERFORMANCE  -DUSE_GRACKLE  -DUSE_UUID -DSAB\n");
   fprintf (fp,"\n");
   fprintf (fp,"INCLUDES = -I/home/jhl1862/library/hdf5/include  -I/home/jhl1862/library/openmpi-4/include       -I/home/jhl1862/code/grackle-3/include    -I.\n");
   fprintf (fp,"\n");
   fprintf (fp,"CPPFLAGS = -P -traditional \n");
   fprintf (fp,"CFLAGS   =  -O2\n");
   fprintf (fp,"CXXFLAGS =  -O2\n");
   fprintf (fp,"FFLAGS   = -fno-second-underscore -ffixed-line-length-132 -O2\n");
   fprintf (fp,"F90FLAGS = -fno-second-underscore -O2\n");
   fprintf (fp,"LDFLAGS  =  -O2\n");
   fprintf (fp,"\n");
   fprintf (fp,"LIBS     = -L/home/jhl1862/library/hdf5/lib -lhdf5 -lz -lgfortran  -L/home/jhl1862/library/openmpi-4/lib        -L/home/jhl1862/code/grackle-3/lib -lgrackle\n");
   fprintf (fp,"\n");
}
