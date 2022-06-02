/*
	 velocity_lib.c (c)
	 
	 Purpose: Functions to update velocity model.
	 	 
	 Version 1.0
	 
	 Site: https://dirack.github.io
	 
	 Programmer: Rodolfo A. C. Neves (Dirack) 19/09/2021

	 Email:  rodolfo_profissional@hotmail.com

	 License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.

*/

#include <stdio.h>
#include <stdlib.h>
#include <rsf.h>
#include "velocity_lib.h"

void updateVelocityModel(  int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velocity model disturbance */
			   int nsv, /* Dimension of sv the vector */
			   float *vel, /* Velocity model */
			   int nvel /* Dimension of the vel vector */)
/*< Velocity model update
Note: This function uses a sv (layers velocity) vector and sz (depth interfaces
coordinates) vector to build the depth velocity model. There is nsv constant
velocity layers in the model and nsv-1 interfaces separating them.
These interfaces are described with nsz control points (nodes) in the sz vector and
they are interpolated using natural cubic spline interpolation.
 >*/
{
        sf_eno2 map;
        float f2[2];
        int i, j, i1, i2; 
        float x, y;
	float *vprof;

	vprof = sf_floatalloc(nsv*4);

	for(j=0;j<4;j++){
		for(i=0;i<nsv;i++){
			vprof[i+nsv*j]=sv[i];
		}
	}

        map = sf_eno2_init(3,nsv,4);

        sf_eno2_set1(map,vprof);
        for(i2=0;i2<n[1];i2++){

                for(i1=0;i1<n[0];i1++){
                        x = (i1*d[0]+o[0])/0.5; i=x; x -= i;
                        y = (i2*d[1]+o[1])/0.5; j=y; y -= j;
                        sf_eno2_apply(map,i,j,x,y,&vel[i2*n[0]+i1],f2,FUNC);
                }
        }
        sf_eno2_close(map);
	free(vprof);
}

void buildSlownessModelFromVelocityModel(int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			 		 float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
					 float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
					 float *sv, /* Velociy disturbance */
					 int nsv, /* Dimension of sv vector */
					 float *vel, /* Velocity model */
					 int nslow /* Dimension of vel vector */)
/*< Slowness model build from velocity model
Note: This function is a function wrapper to updateVelocityModel function.
It calls that function to update the velocity model and build the slowness
model matrix using the slowness definition slow=(1.0/(v*v)). 
 >*/
{

	int i, nm; // Loop counters and indexes

	nm =n[0]*n[1];
	updateVelocityModel(n,o,d,sv,nsv,vel,nm);

	/* transform velocity to slowness */
	for(i=0;i<nm;i++){
			vel[i] = 1.0/(vel[i]*vel[i]);
	}
}

