#include <pthread.h>
#include <stdio.h>
#include <assert.h>
#include "thresfilter.h"

#define NUM_THREADS 5

#define uint unsigned int 

typedef struct {
  int thread_id;
  int xsize;
  int ysize;
  pixel* src;
  int sum;
} thread_args;

void *calculate_average(void *argument) {
  thread_args* ta = (thread_args*)argument;

  int nump = ta->xsize * ta->ysize;
  int i;
  for(i = 0, ta->sum = 0; i < nump; i++) {
    ta->sum += (uint)ta->src[i].r + (uint)ta->src[i].g + (uint)ta->src[i].b;
  }

  printf("%d: %d\n", ta->thread_id, ta->sum);

  return NULL;
}

void thresfilter(const int xsize, const int ysize, pixel* src){


  uint sum, i, psum, nump;

  nump = xsize * ysize;

  pthread_t thread[NUM_THREADS];
  thread_args thread_args[NUM_THREADS];

  int rows = ysize / NUM_THREADS + 1;

  int rc;

  for(i = 0; i < NUM_THREADS; ++i) {
    thread_args[i].thread_id = i;
    thread_args[i].src = src + i * rows * xsize;
    thread_args[i].xsize = xsize;
    thread_args[i].ysize = rows;
    thread_args[i].sum = 0;

    rc = pthread_create(&thread[i], NULL, calculate_average, (void*)&thread_args[i]);
    assert(0 == rc);
  }

  int totsum = 0;
  for(i = 0; i < NUM_THREADS; ++i) {
    rc = pthread_join(thread[i], NULL);
    assert(0 == rc);
    totsum += thread_args[i].sum;
  }

  printf("totsum: %d\n", totsum);
  sum = totsum / nump;
  /*
  for(i = 0, sum = 0; i < nump; i++) {
    sum += (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
  }

  sum /= nump;
  */
  for(i = 0; i < nump; i++) {
    psum = (uint)src[i].r + (uint)src[i].g + (uint)src[i].b;
    if(sum > psum) {
      src[i].r = src[i].g = src[i].b = 0;
    }
    else {
      src[i].r = src[i].g = src[i].b = 255;
    }
    }
}
