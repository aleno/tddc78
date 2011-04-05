#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include <mpi.h>
#include "ppmio.h"
#include "blurfilter.h"
#include "gaussw.h"


#define ROOT 0

int range(int value, int min, int max)
{
  if(value < min)
    return min;
  else if(value > max)
    return max;
  return value;
}

int main (int argc, char ** argv) {
   int radius;
    int xsize, ysize, colmax;
    pixel src[MAX_PIXELS];
    pixel dst[MAX_PIXELS];
    struct timespec stime, etime;
#define MAX_RAD 1000

    double w[MAX_RAD];

    /* Take care of the arguments */

    if (argc != 4) {
	fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
	exit(1);
    }
    radius = atoi(argv[1]);
    if((radius > MAX_RAD) || (radius < 1)) {
	fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
	exit(1);
    }
    int numprocs;
    int myid;
    int i;
 
    MPI_Init(&argc,&argv); /* all MPI programs start with MPI_Init; all 'N' processes exist thereafter */
    MPI_Comm_size(MPI_COMM_WORLD,&numprocs); /* find out how big the SPMD world is */
    MPI_Comm_rank(MPI_COMM_WORLD,&myid); /* and this processes' rank is */
    
    if( myid == ROOT) {
      /* read file */
      if(read_ppm (argv[2], &xsize, &ysize, &colmax, (char *) src) != 0)
        exit(1);
      
      if (colmax > 255) {
	fprintf(stderr, "Too large maximum color-component value\n");
	exit(1);
      }
    }


    int sizes[2];
    sizes[0] = xsize;
    sizes[1] = ysize;

    MPI_Bcast(sizes, 2, MPI_INT, ROOT, MPI_COMM_WORLD);

    if( myid != ROOT) {
      xsize = sizes[0];
      ysize = sizes[1];
    }

    // Skapa buffrar för gatherv/scatterv
    int* displs = (int*)malloc(sizeof(int) * numprocs);
    int* scounts = (int*)malloc(sizeof(int) * numprocs);
    int* srows = (int*)malloc(sizeof(int) * numprocs);

    int rows = (ysize / numprocs);
    
    // Största displacement för bilden.
    int max_size = ysize * xsize * sizeof(pixel);

    for(i = 0; i < numprocs; i++) {
      displs[i] = (i * rows - radius) * xsize * sizeof(pixel);
      srows[i] = (2 * radius + rows); // * xsize * sizeof(pixel);
      scounts[i] = srows[i] * xsize * sizeof(pixel);

      // Fixa det dumma fallet med radius > rows per node.
      //displs[i] = range(displs[i], 0, max_size);
      //scounts[i] = range(displs[i] + scounts[i], 0, max_size) - displs[i];

      /*      if(displs[i] + scounts[i] > ysize * xsize * sizeof(pixel))
	scounts[i] = ysize* xsize * sizeof(pixel) - displs[i];
      if(displs[i] < 0) {
	scounts[i] = ysize* xsize * sizeof(pixel) - displs[i];
	displs[i] = 0;
	}*/
    }

    // Specialfall för första/sista noden
    displs[0] = 0;
    srows[0] = radius + rows;
    srows[numprocs - 1] = (ysize - (numprocs - 1) * rows) + radius;
    scounts[0] = srows[0] * xsize * sizeof(pixel);
    scounts[numprocs - 1] = (ysize - srows[numprocs - 1]) * xsize * sizeof(pixel);

    MPI_Scatterv(src, scounts, displs, MPI_BYTE,
		 dst, scounts[myid], MPI_BYTE,
		 ROOT, MPI_COMM_WORLD);

    //    if(myid == ROOT)
    //  memcpy(dst, src, sizeof(pixel) * MAX_PIXELS);

    printf("Task %d, Has read the image, generating coefficients\n", myid);

    /* filter */
    get_gauss_weights(radius, w);

    printf("Task %d, Calling filter\n", myid);

    clock_gettime(CLOCK_REALTIME, &stime);

    printf("Task %d, X %d, Y %d, R %d\n", myid, xsize, ysize, srows[myid]);

    blurfilter(xsize, srows[myid], dst, radius, w);

    printf("Task %d, Has run filter, Gather next\n", myid);
        clock_gettime(CLOCK_REALTIME, &etime);

    printf("Task %d: Filtering took: %g secs\n", myid, (etime.tv_sec  - stime.tv_sec) +
       1e-9*(etime.tv_nsec  - stime.tv_nsec)) ;

    pixel* buffp = dst + radius * xsize;
    if(myid == ROOT)
      buffp = dst;

    for(int i = 0; i < numprocs; i++) {
      srows[i] = rows * xsize * sizeof(pixel);
      displs[i] = rows * i * xsize * sizeof(pixel);
      scounts[i] = rows * xsize * sizeof(pixel);
    }
    scounts[numprocs - 1] = (ysize - rows * (numprocs - 1)) * xsize * sizeof(pixel);

    MPI_Gatherv(buffp, scounts[myid], MPI_BYTE,
		src, scounts, displs, MPI_BYTE,
	       ROOT, MPI_COMM_WORLD);
    printf("Task %d, Has gather data\n", myid);

    if(myid == ROOT) {
      /* write result */
      printf("Task %d, Writing output file\n", myid);
      
      if(write_ppm (argv[3], xsize, ysize, (char *)src) != 0) {
	printf("Didnt work!!\n");
	exit(1);
      }
      printf("Written!!!!\n");
    }

    printf("Task %d, Terminate\n", myid);

    MPI_Finalize(); /* MPI Programs end with MPI Finalize; this is a weak synchronization point */



    return(0);
}
