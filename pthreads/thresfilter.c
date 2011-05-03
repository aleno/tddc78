#include <pthread.h>
#include <semaphore.h>
#include <stdio.h>
#include <assert.h>
#include "thresfilter.h"

#define NUM_THREADS 8

#define uint unsigned int 

typedef struct {
  int thread_id;
  int xsize;
  int ysize;
  int rows;
  pixel* src;
} thread_args;

pthread_mutex_t stage_lock;
pthread_cond_t stage_cond;

int globsum;
int threads_in_stage = 0;

void *calculate_average(void *argument) {
  thread_args* ta = (thread_args*)argument;

  int nump = ta->xsize * ta->rows;
  int i, sum;
  for(i = 0, sum = 0; i < nump; i++) {
    sum += (uint)ta->src[i].r + (uint)ta->src[i].g + (uint)ta->src[i].b;
  }

  printf("%d: %d\n", ta->thread_id, sum);
  pthread_mutex_lock(&stage_lock);
  globsum += sum;
  threads_in_stage --;
  while( threads_in_stage > 0 ) 
    pthread_cond_wait(&stage_cond, &stage_lock);
  pthread_cond_signal(&stage_cond);
  pthread_mutex_unlock(&stage_lock);

  printf("%d: globsum %d\n", ta->thread_id, globsum);
  
  sum = globsum / (ta->xsize * ta->ysize);
  int psum;
  for(i = 0; i < nump; i++) {
    psum = (uint)ta->src[i].r + (uint)ta->src[i].g + (uint)ta->src[i].b;
    if(sum > psum) {
      ta->src[i].r = ta->src[i].g = ta->src[i].b = 0;
    }
    else {
      ta->src[i].r = ta->src[i].g = ta->src[i].b = 255;
    }
  }

  return NULL;
}

void thresfilter(const int xsize, const int ysize, pixel* src){


  uint sum, i, psum, nump;

  nump = xsize * ysize;

  pthread_t thread[NUM_THREADS];
  thread_args thread_args[NUM_THREADS];

  int rows = ysize / NUM_THREADS + 1;

  threads_in_stage = NUM_THREADS;

  int rc;

  for(i = 0; i < NUM_THREADS; ++i) {
    thread_args[i].thread_id = i;
    thread_args[i].src = src + i * rows * xsize;
    thread_args[i].xsize = xsize;
    thread_args[i].ysize = ysize;
    thread_args[i].rows = rows;

    rc = pthread_create(&thread[i], NULL, calculate_average, (void*)&thread_args[i]);
    assert(0 == rc);
  }

  for(i = 0; i < NUM_THREADS; ++i) {
    rc = pthread_join(thread[i], NULL);
    assert(0 == rc);
  }
}
