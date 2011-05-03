/*
  File: blurfilter.c

  Implementation of blurfilter function.
    
 */
#include <stdio.h>
#include "blurfilter.h"
#include "ppmio.h"

#include <pthread.h>
#include <semaphore.h>
#include <assert.h>

extern int NUM_THREADS;

pixel* pix(pixel* image, const int xx, const int yy, const int xsize)
{
  register int off = xsize*yy + xx;

#ifdef DBG
  if(off >= MAX_PIXELS) {
    fprintf(stderr, "\n Terribly wrong: %d %d %d\n",xx,yy,xsize);
  }
#endif
  return (image + off);
}

struct blurfilter_args {
  int thread_id;
  int xsize;
  int ysize;
  int rows;
  int radius;
  const double *w;
  pixel* src;
};

int reader_threads;

pthread_mutex_t reader_lock;
pthread_cond_t reader_cond;

#define min(a,b) ((a<b)?a:b)
#define max(a,b) ((a>b)?a:b)

void *blurfilter_thread(void *args) {
  struct blurfilter_args *ta = (struct blurfilter_args*)args;

  int x,y,x2,y2, wi;
  double r,g,b,n, wc;
  pixel dst[MAX_PIXELS];

  const int xsize = ta->xsize;
  const int ysize = ta->ysize;
  pixel *src = ta->src;
  const double *w = ta->w;
  const int radius = ta->radius;

  int start = max(0, ta->thread_id * (ta->rows - radius));
  int end = min((ta->thread_id + 1) * (ta->rows + radius), ysize);

  //  printf("%d: do %d..%d\n", ta->thread_id, start, end);

  for (y=start; y<end; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(src, x, y, xsize)->r;
      g = w[0] * pix(src, x, y, xsize)->g;
      b = w[0] * pix(src, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
	wc = w[wi];
	x2 = x - wi;
	if(x2 >= 0) {
	  r += wc * pix(src, x2, y, xsize)->r;
	  g += wc * pix(src, x2, y, xsize)->g;
	  b += wc * pix(src, x2, y, xsize)->b;
	  n += wc;
	}
	x2 = x + wi;
	if(x2 < xsize) {
	  r += wc * pix(src, x2, y, xsize)->r;
	  g += wc * pix(src, x2, y, xsize)->g;
	  b += wc * pix(src, x2, y, xsize)->b;
	  n += wc;
	}
      }
      pix(dst,x,y, xsize)->r = r/n;
      pix(dst,x,y, xsize)->g = g/n;
      pix(dst,x,y, xsize)->b = b/n;
    }
  }

  // Reader writer problem..
  // Wait for all readers to finish before letings writers in..
  
  //  printf("%d: Hold and wait...", ta->thread_id);
  pthread_mutex_lock(&reader_lock);
  reader_threads --;
  while( reader_threads > 0 ) 
    pthread_cond_wait(&reader_cond, &reader_lock);
  pthread_cond_signal(&reader_cond);
  pthread_mutex_unlock(&reader_lock);

  //  printf("%d: Gooo...\n", ta->thread_id);
  start = ta->thread_id * ta->rows;
  end = min((ta->thread_id + 1) * (ta->rows), ysize);

  for (y=start; y<end; y++) {
    for (x=0; x<xsize; x++) {
      r = w[0] * pix(dst, x, y, xsize)->r;
      g = w[0] * pix(dst, x, y, xsize)->g;
      b = w[0] * pix(dst, x, y, xsize)->b;
      n = w[0];
      for ( wi=1; wi <= radius; wi++) {
	wc = w[wi];
	y2 = y - wi;
	if(y2 >= 0) {
	  r += wc * pix(dst, x, y2, xsize)->r;
	  g += wc * pix(dst, x, y2, xsize)->g;
	  b += wc * pix(dst, x, y2, xsize)->b;
	  n += wc;
	}
	y2 = y + wi;
	if(y2 < ysize) {
	  r += wc * pix(dst, x, y2, xsize)->r;
	  g += wc * pix(dst, x, y2, xsize)->g;
	  b += wc * pix(dst, x, y2, xsize)->b;
	  n += wc;
	}
      }
      pix(src,x,y, xsize)->r = r/n;
      pix(src,x,y, xsize)->g = g/n;
      pix(src,x,y, xsize)->b = b/n;
    }
  }

  return NULL;
}

void blurfilter(const int xsize, const int ysize, pixel* src, const int radius, const double *w){
  
  pthread_t thread[NUM_THREADS];
  struct blurfilter_args thread_args[NUM_THREADS];

  int rows = ysize / NUM_THREADS + 1;

  reader_threads = NUM_THREADS;

  int rc, i;

  printf("Start %d threads...\n", reader_threads);

  for(i = 0; i < NUM_THREADS; ++i) {
    thread_args[i].thread_id = i;
    thread_args[i].src = src;
    thread_args[i].xsize = xsize;
    thread_args[i].ysize = ysize;
    thread_args[i].rows = rows;
    thread_args[i].radius = radius;
    thread_args[i].w = w;

    rc = pthread_create(&thread[i], NULL, blurfilter_thread,
			(void*)&thread_args[i]);
    assert(0 == rc);
  }

  for(i = 0; i < NUM_THREADS; ++i) {
    rc = pthread_join(thread[i], NULL);
    assert(0 == rc);
  }
}
