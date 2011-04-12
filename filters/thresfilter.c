#include <mpi.h>
#include <stdio.h>
#include <malloc.h>
#include "thresfilter.h"

extern int myid;
extern int numprocs;

void thresfilter(const int xsize, const int ysize, pixel* src){
#define uint unsigned int 

  uint sum, i, psum, nump;

  nump = xsize * ysize;

  int rows = ysize / numprocs;

  pixel dst[1000000];

  MPI_Scatter(src, rows * xsize * sizeof(pixel), MPI_BYTE,
	      dst, rows * xsize * sizeof(pixel), MPI_BYTE,
	      0, MPI_COMM_WORLD);

  int start = myid * rows * xsize;
  int end = (myid + 1) * rows * xsize;
  if(end > nump) end = nump;

  //for(i = start, sum = 0; i < end; i++) {
  for(i = 0; i < rows * xsize; i++) {
    sum += (uint)dst[i].r + (uint)dst[i].g + (uint)dst[i].b;
  }

  printf("Task %d, Sum is %d\n", myid, sum);

  //sum /= end - start;
  //int* sum_v = (int*)malloc(sizeof(int) * numprocs);

  int totsum = 0;

  MPI_Reduce( &sum, &totsum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

  //totsum /= ysize * xsize;
  sum = totsum / (ysize * xsize);

  MPI_Bcast(&totsum, 1, MPI_INT, 0, MPI_COMM_WORLD);

  for(i = 0; i < nump; i++) {
//  for(i = start; i < end; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(sum > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
  }
}
