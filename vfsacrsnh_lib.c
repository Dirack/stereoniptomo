/*
	 vfsacrsnh_lib.c (c)
	 
	 Purpose: 'Mvfsacrsnh.c' library.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2019

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

/*
TODO: Modify macro definition in search window for each interface.
Large windows can make the result oscilate a lot and do not converge
*/
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>
/*^*/

#define signal(s) ((s<0)?(-1.):(1.))
/*< Signal function >*/
/*^*/

float getRandomNumberBetween0and1(){
/*< Function to get a random number between 0 and 1 >*/

	return (float)(rand()%1000)/1000;
}

float getVfsaIterationTemperature(int iteration,float dampingFactor,float inicialTemperature){
/*< Temperature function for VFSA algorithm >*/

	return inicialTemperature*expf(-dampingFactor*pow(iteration,0.25));

}

/* TODO: Modify this function for multiple interfaces */
void disturbParameters( float temperature, /* Temperature of this interation in VFSA */
			float* disturbedVel, /* Parameters disturbed vector */
			float *vel,
			float minvel,
			float maxvel,
			float scale /* Scale to multiply by disturbance */)
/*< Disturb parameters from the previous iteration of VFSA
 Note: It receives a parameter vector and distubs it accordingly to 
VFSA disturb parameters step.
 >*/
{

	float u;
	float disturbance;
	int i;
	float min[6]={1.5,1.5,1.65,1.65,1.75,2.0};
	float max[6]={1.5,1.65,1.75,1.85,2.0,3.0};

	disturbedVel[0]=vel[0];
	disturbedVel[5]=vel[5];
	for(i=1;i<5;i++){

		u=getRandomNumberBetween0and1();
					
		disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);
		disturbance *= scale;

		disturbedVel[i] = vel[i] + (disturbance);

		if (disturbedVel[i] >= max[i])
			disturbedVel[i] = max[i] - (max[i]-min[i]) * getRandomNumberBetween0and1();

		if (disturbedVel[i] <= min[i])
			disturbedVel[i] = (max[i]-min[i]) * getRandomNumberBetween0and1() + min[i];
	}
}

/* TODO: Modify this function for multiple interfaces */
void disturbParameters2( float temperature, /* Temperature of this interation in VFSA */
			float *disturbedVel, /* Parameters disturbed vector */
			float *vel,
			int *nsv,
			float minvel,
			float maxvel,
			float scale /* Scale to multiply by disturbance */)
/*< Disturb parameters from the previous iteration of VFSA
 Note: It receives a parameter vector and distubs it accordingly to 
VFSA disturb parameters step.
 >*/
{

	float u;
	float disturbance;
	int i,j;
	float min[6]={1.5,1.5,1.65,1.65,1.75,2.0};
	float max[6]={1.5,1.65,1.75,1.85,2.0,3.0};

	for(j=0;j<nsv[1];j++){
		disturbedVel[0+j*nsv[0]]=vel[0+j*nsv[0]];
		disturbedVel[nsv[0]-1+j*nsv[0]]=vel[nsv[0]-1+j*nsv[0]];
	}

	for(j=0;j<nsv[1];j++){
		for(i=1;i<nsv[0]-1;i++){

			u=getRandomNumberBetween0and1();
						
			disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);
			disturbance *= scale;

			disturbedVel[i+j*nsv[0]] = vel[i+j*nsv[0]] + (disturbance);

			if (disturbedVel[i+j*nsv[0]] >= max[i])
				disturbedVel[i+j*nsv[0]] = max[i] - (max[i]-min[i]) * getRandomNumberBetween0and1();

			if (disturbedVel[i+j*nsv[0]] <= min[i])
				disturbedVel[i+j*nsv[0]] = (max[i]-min[i]) * getRandomNumberBetween0and1() + min[i];
		}
	}
	//dumpfloat1("dis",disturbedVel,nsv[0]*nsv[1]);
	//sf_error("oi");
}

/* TODO: Modify this function for multiple interfaces */
void disturbParameters3( float temperature, /* Temperature of this interation in VFSA */
			float *disturbedVel, /* Parameters disturbed vector */
			float *vel,
			float minvel,
			float maxvel,
			float scale /* Scale to multiply by disturbance */)
/*< Disturb parameters from the previous iteration of VFSA
 Note: It receives a parameter vector and distubs it accordingly to 
VFSA disturb parameters step.
 >*/
{

	float u;
	float disturbance;
	int i,j;
	float min[6]={1.5,1.5,1.65,1.65,1.75,2.0};
	float max[6]={1.5,1.65,1.75,1.85,2.0,3.0};

	//for(j=0;j<nsv[1];j++){
	//	disturbedVel[0+j*nsv[0]]=vel[0+j*nsv[0]];
	//	disturbedVel[nsv[0]-1+j*nsv[0]]=vel[nsv[0]-1+j*nsv[0]];
	//}

	//for(j=0;j<nsv[1];j++){
	//	for(i=1;i<nsv[0]-1;i++){

			u=getRandomNumberBetween0and1();
						
			disturbance = signal(u - 0.5) * temperature * (pow( (1+temperature),fabs(2*u-1) )-1);
			disturbance *= scale;

			disturbedVel[0] = vel[0] + (disturbance);

			if (disturbedVel[0] >= maxvel)
				disturbedVel[0] = maxvel - (maxvel-minvel) * getRandomNumberBetween0and1();

			if (disturbedVel[0] <= minvel)
				disturbedVel[0] = (maxvel-minvel) * getRandomNumberBetween0and1() + minvel;
	//	}
	//}
	//dumpfloat1("dis",disturbedVel,nsv[0]*nsv[1]);
	//sf_error("oi");
}

