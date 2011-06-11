#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include <mpi.h>
#include "definitions.h"

#define ROOT 0

int main (int argc, char ** argv) {
  /*if (argc != 4) {
    fprintf(stderr, "Usage: %s radius infile outfile\n", argv[0]);
    exit(1);
    }
    radius = atoi(argv[1]);
    if((radius > MAX_RAD) || (radius < 1)) {
    fprintf(stderr, "Radius (%d) must be greater than zero and less then %d\n", radius, MAX_RAD);
    exit(1);
    }*/
  
  int numprocs;
  int myid;
  int i;
 
  MPI_Init(&argc,&argv); /* all MPI programs start with MPI_Init; all 'N' processes exist thereafter */
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs); /* find out how big the SPMD world is */
  MPI_Comm_rank(MPI_COMM_WORLD,&myid); /* and this processes' rank is */

  int sub_box_width=BOX_HORIZ_SIZE/numprocs;
  int sub_box_height=BOX_VERT_SIZE;

  cord_t big_box;
  big_box.x0 = 0;
  big_box.x1 = BOX_HORIZ_SIZE;
  big_box.y0 = 0;
  big_box.y1 = BOX_VERT_SIZE;

  //Calculate the box of this processor
  cord_t box;
  box.x0 = myid*sub_box_width;
  box.x1 = (myid+1)*sub_box_width;
  box.y0 = myid*sub_box_height;
  box.y1 = (myid+1)*sub_box_height;

  //Initiate particles
  int numparticles = atoi(argv[1]);
  int active_particles = numparticles / numprocs;
  numparticles = active_particles * numprocs;

  struct particle *particles =
    (struct particle*)malloc(sizeof(struct particle) * numparticles);
  struct particle *to_east =
    (struct particle*)malloc(sizeof(struct particle) * numparticles);
  int east_emigrants = 0;
  struct particle *to_west =
    (struct particle*)malloc(sizeof(struct particle) * numparticles);
  int west_emigrants = 0;

  for (int i = 0; i < active_particles; i++) {
    int x = rand() % sub_box_width;
    int y = rand() % sub_box_height;
      
    particles[i].pcord.x = box.x0 + x;
    particles[i].pcord.y = box.y0 + y;
    int vel = rand() % 50;
    int dir = (rand() % 360) * 2.0/360.0*PI;
    particles[i].pcord.vx = vel * cos(dir);
    particles[i].pcord.vy = vel * sin(dir);
  }     
      

  if( myid == ROOT) {
      
  }

  float time_step = 0.5;
  float momentum = 0;
  int steps = 100;

  //Main loop: for each time-step do
  while (steps-- > 0) {
    printf("Steps left: %d\n", steps);
      
    //For all particles do:
    for(int i = 0; i < active_particles; i++) {
      int collided = 0;
      for(int j = i + 1; j < active_particles; j++) {
	//Check for collisions
	float diff = collide(&particles[i].pcord, &particles[j].pcord);
	if(diff != -1) {
	  collided = 1;
	  interact(&particles[i].pcord, &particles[j].pcord, diff);
	}
      }
      if(collided == 0) {
	//Move particles that has not collided with another
	feuler(&particles[i].pcord, time_step);
      }

      //Chech for wall interaction and add the momentum
      momentum += wall_collide(&particles[i].pcord, big_box);

      // Add particle to move list of out of box. =/
      if(particles[i].pcord.x >= box.x1){
	to_east[east_emigrants++] = particles[i];
	particles[i] = particles[--active_particles];
      } else if(particles[i].pcord.x < box.x0){
	to_west[west_emigrants++] = particles[i];
	particles[i] = particles[--active_particles];
	// Swap with last.
	// active --
	// to_move ++;
      }
      // Do that..
    }

    //Communicte particles Aouwvhh.      
    // Notify goverments about new imigrants.. =/
    /*int to_send = east_emigrants;
        MPI_Send(&to_send, 1, MPI_INT, (myid + 1) % numprocs,
	     myid, 99, MPI_COMM_WORLD);
	
    int to_recv;
    MPI_Recv(&to_recv, 1, MPI_INT, (myid + numprocs - 1) % numprocs, 99, MPI_COMM_WORLD);*/
		
    int local[2] = {east_emigrants, west_emigrants};
    int ls[2] = {0};

    int lnbr = (myid + numprocs - 1) % numprocs;
    int rnbr = (myid + 1) % numprocs;

    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Status status;

    // Exchange boundary values with neighbors:
    MPI_Send( local+1, 1, MPI_INT, lnbr, 10, com );
    MPI_Recv( ls+1, 1, MPI_INT, rnbr, 10, com, &status );
    MPI_Send( ls, 1, MPI_INT, rnbr, 20, com );
    MPI_Recv( local, 1, MPI_INT, lnbr, 20, com, &status );
    // Let them flee...
      
    MPI_Send( to_west, local[1]* sizeof(struct particle), MPI_BYTE,
	      lnbr, 30, com );
    MPI_Recv( particles+active_particles, ls[1] * sizeof(struct particle),
	      MPI_BYTE, rnbr, 30, com, &status );
    active_particles += ls[1];
      
    MPI_Send( to_east, local[1]* sizeof(struct particle), MPI_BYTE,
	      rnbr, 40, com );
    MPI_Recv( particles+active_particles, ls[0] * sizeof(struct particle),
	      MPI_BYTE, lnbr, 40, com, &status );
    active_particles += ls[0];

    // Reduce momentum
    float momentum2;
    MPI_Allreduce(&momentum, &momentum2, 1, MPI_FLOAT, MPI_SUM, com);
    momentum = momentum2;
    // Show what to do..
    printf("Momentum: %f\n", momentum);

    // Maybe say that to all..
  }

  // Fin, :D
    
  //Calculate pressure
  /*
  double starttime = MPI_Wtime();


  MPI_Bcast(sizes, 2, MPI_INT, ROOT, MPI_COMM_WORLD);

  if( myid != ROOT) {
      
  }

  // Skapa buffrar för gatherv/scatterv
  int* displs = (int*)malloc(sizeof(int) * numprocs);
  int* scounts = (int*)malloc(sizeof(int) * numprocs);
  int* srows = (int*)malloc(sizeof(int) * numprocs);

  int rows = (ysize / numprocs);
    
  for(i = 0; i < numprocs; i++) {

  }

  MPI_Scatterv(src, scounts, displs, MPI_BYTE,
	       dst, scounts[myid], MPI_BYTE,
	       ROOT, MPI_COMM_WORLD);

  MPI_Gatherv(buffp, scounts[myid], MPI_BYTE,
	      src, scounts, displs, MPI_BYTE,
	      ROOT, MPI_COMM_WORLD);
    */
  MPI_Finalize(); /* MPI Programs end with MPI Finalize; this is a weak synchronization point */

  return(0);
}

