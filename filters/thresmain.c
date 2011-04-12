#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include "ppmio.h"
#include "thresfilter.h"

int numprocs;
int myid;

int main (int argc, char ** argv) {
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    struct timespec stime, etime;

    /* Take care of the arguments */

    if (argc != 3) {
	fprintf(stderr, "Usage: %s infile outfile\n", argv[0]);
	exit(1);
    }

    MPI_Init(&argc,&argv); /* all MPI programs start with MPI_Init; all 'N' processes exist thereafter */
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); /* find out how big the SPMD world is */
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); /* and this processes' rank is */
    
    if(myid == 0) {
      /* read file */
      if(read_ppm (argv[1], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);
      
      if (colmax > 255) {
	fprintf(stderr, "Too large maximum color-component value\n");
	exit(1);
      }
    }

    int sizes[2];
    sizes[0] = xsize;
    sizes[1] = ysize;

    MPI_Bcast(sizes, 2, MPI_INT, 0, MPI_COMM_WORLD);

    if( myid != 0) {
      xsize = sizes[0];
      ysize = sizes[1];
    }

    printf("Task %d, Size is %d x %d\n", myid, xsize, ysize);


    printf("Has read the image, calling filter\n");

    //clock_gettime(CLOCK_REALTIME, &stime);

    thresfilter(xsize, ysize, src);

    //clock_gettime(CLOCK_REALTIME, &etime);

    // printf("Filtering took: %g secs\n", (etime.tv_sec  - stime.tv_sec) +
    //   1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    /* write result */
    if(myid == 0) {
      printf("Writing output file\n");
      
      if(write_ppm (argv[2], xsize, ysize, (char *)src) != 0)
	exit(1);
    }

    MPI_Finalize(); /* MPI Programs end with MPI Finalize; this is a weak synchronization point */



    return(0);
}
