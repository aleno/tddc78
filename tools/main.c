#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <malloc.h>
#include <time.h>
#include <mpi.h>
#include <VT.h>
#include "definitions.h"

#define ROOT 0
#define NUM_STEPS 100

int application_class;

int main (int argc, char ** argv) {
  int numprocs;
  int myid;
  int i;
 
  srand(time(0));

  /* all MPI programs start with MPI_Init;
     all 'N' processes exist thereafter */
  MPI_Init(&argc,&argv);
  /* find out how big the SPMD world is */
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  /* and this processes' rank is */
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  VT_classdef ("Partsim Application", &application_class);

  if(myid == ROOT) {
    printf("Number of processors: %d\n", numprocs);
  }

  int migration_state;
  int reduce_state;
  int step_state;

  VT_funcdef("Migration", application_class, &migration_state);
  VT_funcdef("Reduce", application_class, &reduce_state);
  VT_funcdef("Step", application_class, &step_state);

  int pd_counter_handle;
  int counter_class;
  const float boundaries[2]={0,atoi(argv[1])};
  
  VT_countdef( "Particles", application_class,
	       VT_COUNT_INTEGER|VT_COUNT_ABSVAL|VT_COUNT_VALID_AFTER,
	       VT_ME,
	       boundaries,"#",&pd_counter_handle);
  
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

  if(myid == ROOT)
    printf("There shall be %d particles.\n", numparticles);

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
  int steps = NUM_STEPS;
  double start_time = MPI_Wtime();

  //Main loop: for each time-step do
  while (steps-- > 0) {
    VT_countval(1, &pd_counter_handle, &active_particles);


    VT_enter(step_state, VT_NOSCL);
    //printf("Steps left: %d\n", steps);
      
    float momentum_tmp = 0;
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
      momentum_tmp += wall_collide(&particles[i].pcord, big_box);

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
    VT_leave(VT_NOSCL);

    VT_enter(migration_state, VT_NOSCL);
    //Communicte particles.
    int local[2] = {east_emigrants, west_emigrants};
    int ls[2] = {0};

    int lnbr = (myid + numprocs - 1) % numprocs;
    int rnbr = (myid + 1) % numprocs;

    MPI_Comm com = MPI_COMM_WORLD;
    MPI_Status status;

    // Exchange boundary values with neighbors:
    MPI_Send( local+1, 1, MPI_INT, lnbr, 10, com );
    MPI_Recv( ls+1, 1, MPI_INT, rnbr, 10, com, &status );
    MPI_Send( local, 1, MPI_INT, rnbr, 20, com );
    MPI_Recv( ls, 1, MPI_INT, lnbr, 20, com, &status );
    // Let them flee...
      
    MPI_Send( to_west, local[1]* sizeof(struct particle), MPI_BYTE,
	      lnbr, 30, com );
    MPI_Recv( particles+active_particles, ls[1] * sizeof(struct particle),
	      MPI_BYTE, rnbr, 30, com, &status );
    active_particles += ls[1];
      
    MPI_Send( to_east, local[0]* sizeof(struct particle), MPI_BYTE,
	      rnbr, 40, com );
    MPI_Recv( particles+active_particles, ls[0] * sizeof(struct particle),
	      MPI_BYTE, lnbr, 40, com, &status );
    active_particles += ls[0];

    if(local[0] != 0 || local[1] != 0 || ls[0] != 0 || ls[1] != 0)
    printf("%d: %d/%d - %d/%d %d\n", myid,
	   local[0], local[1], ls[0], ls[1], active_particles);


    west_emigrants = 0;
    east_emigrants = 0;

    VT_leave( VT_NOSCL );

    momentum += momentum_tmp;

    // Maybe say that to all..
  }
  
  printf("%d: Pre momentum: %f, %d\n", myid, momentum, active_particles);
  VT_enter(reduce_state, VT_NOSCL);
  // Reduce momentum
  float momentum2;
  MPI_Reduce(&momentum, &momentum2, 1, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD);
  momentum = momentum2;
  // Show what to do..
  //printf("Momentum: %f\n", momentum);
  VT_leave(VT_NOSCL);
  double end_time = MPI_Wtime();

  int foo;

  MPI_Reduce(&active_particles, &foo, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

  momentum /= NUM_STEPS * (2*BOX_VERT_SIZE + 2* BOX_HORIZ_SIZE);

  if(myid == ROOT) {
    printf("The mometum ended up as... %f\n", momentum);
    printf("It's done it took: %lf\n", end_time - start_time);

    printf("Particles in system: %d\n", foo);
  }
  
  // Fin, :D
  /* MPI Programs end with MPI Finalize;
     this is a weak synchronization point */

  MPI_Finalize();

  return(0);
}

