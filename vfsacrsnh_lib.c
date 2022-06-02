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

