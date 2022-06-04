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
#include "raytrace.h"

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

void updateVelocityModel2(
			   float *vel,
			   int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velocity model disturbance */
			   int *nsv, /* Dimension of sv the vector */
			   float *osv,
			   float *dsv)
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

	vprof=sf_floatalloc(nsv[0]*nsv[1]);
	for(i=0;i<(nsv[0]*nsv[1]);i++)
		vprof[i]=sv[i];

	//dumpfloat1("velprof",vprof,nsv[0]*nsv[1]);
        map = sf_eno2_init(3,nsv[0],nsv[1]);

        sf_eno2_set1(map,vprof);
        for(i2=0;i2<n[1];i2++){

                for(i1=0;i1<n[0];i1++){
                        x = (i1*d[0]+o[0]-osv[0])/dsv[0]; i=x; x -= i;
                        y = (i2*d[1]+o[1]-osv[1])/dsv[1]; j=y; y -= j;
                        sf_eno2_apply(map,i,j,x,y,&vel[i2*n[0]+i1],f2,FUNC);
                }
        }
        sf_eno2_close(map);
	free(vprof);
	//sf_warning("%d %d %d %d",nsv[0],nsv[1],n[0],n[1]);
	//dumpfloat1("vel",vel,n[1]*n[0]);
	//sf_error("oi");
}

void buildSlownessModelFromVelocityModel2(
					 float *vel,
					 int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			 		 float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
					 float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
					 float *sv, /* Velociy disturbance */
					 int *nsv, /* Dimension of sv vector */
					 float *osv,
					 float *dsv)
/*< Slowness model build from velocity model
Note: This function is a function wrapper to updateVelocityModel function.
It calls that function to update the velocity model and build the slowness
model matrix using the slowness definition slow=(1.0/(v*v)). 
 >*/
{

	int i, nm; // Loop counters and indexes

	nm =n[0]*n[1];
//	dumpfloat1("antes",vel,n[0]);
	updateVelocityModel2(vel,n,o,d,sv,nsv,osv,dsv);

//	dumpfloat1("depois",vel,n[0]);
	/* transform velocity to slowness */
	for(i=0;i<nm;i++){
			vel[i] = 1.0/(vel[i]*vel[i]);
	}
}

void calcInterfacesZcoord(      float *zi, /* Interfaces depth coordinates */
                                int nint, /* Number of interfaces */
                                float xs, /* x coordinate */
                                int si, /* Spline index */
                                float **coef /* Cubic spline coefficients */)
/*< Calculate depth coordinates of the interfaces
 * Note: This function calculates interfaces depth coordinates and stores it
 * in the zi vector.
  >*/
{
        int i; // Loop counter

        for(i=0;i<nint;i++)
                zi[i] = coef[i][si*4+0]*xs*xs*xs+coef[i][si*4+1]*xs*xs+coef[i][si*4+2]*xs+coef[i][si*4+3];
}

//TODO: correct "coefficients" in the function name
void calculateSplineCoeficients(int n, /* Vectors (x,y) dimension */
                                float* x, /* x coordinates */
                                float* y, /* y coordinates */
                                float* coef /* Spline coefficients */)
/*< Function to calculate natural cubic spline coefficients

Note: It Receives n points and two vectors x and y with n dimension.
It returns a coefficients vector with 4 coefficients for each of the
n-1 natural cubic splines, coef[(n-1)*4].

IMPORTANT: The number of points must be equal or major than 3 (n>3)
and x vector must be in crescent order.

>*/
{

        float s2[n]; // Second derivatives matrix
        int i, ip1, ip2, im1, m; // Loop counter
        float hb, ha, deltaa, deltab, t; // temporary variables
        float e[n-2]; // hi's vector
        float dp[n-2]; // main diagonal

        /* Vectors dimension must be major than 3 */
        if(n<3) sf_error("Vectors dimension n must be major than 3\n");

        /* x vector must be in crescent order */
        for(i=1;i<n;i++){
                if(x[i-1]>x[i]) sf_error("Vector x should be in ascending order\n");
        }

        /* Simetric tridiagonal linear system build */
        ha = x[1]-x[0]; deltaa = (y[1]-y[0])/ha; m=n-2;
        for(i=0;i<m;i++){
                ip1 = i+1; ip2 = i+2;
                hb = x[ip2]-x[ip1];
                deltab = (y[ip2]-y[ip1])/hb;
                e[i] = hb; dp[i] = 2*(ha+hb);
                s2[ip1] = 6*(deltab-deltaa);
                ha=hb; deltaa=deltab;
        }

        /* Gauss elimination */
        for(i=1;i<m;i++){
                ip1=i+1; im1=i-1;
                t = e[im1]/dp[im1];
                dp[i] = dp[i]-t*e[im1];
                s2[ip1] = s2[ip1]-t*s2[i];
        }

        /* Retroactive substitutive solution */
        s2[m]=s2[m]/dp[m-1];
        for(i=m-1;i>0;i--){
                ip1=i+1; im1=i-1;
                s2[i]=(s2[i]-e[im1]*s2[ip1])/dp[im1];
        }
        s2[0]=0; s2[n-1]=0;
        /* Calculate spline coefficients */
        for(i=0;i<n-1;i++){
                ha = x[i+1]-x[i];
                coef[0+i*4] = (s2[i+1]-s2[i])/(6*ha);
                coef[1+i*4] = s2[i]/2;
                coef[2+i*4] = (y[i+1]-y[i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
                coef[3+i*4] = y[i];
        }
}

void interfaceInterpolationFromNipSources(float **s, /* NIP sources */
                                          int ns, /* Number of NIP sources */
                                          float *sz, /* Spline nodepoints */
                                          int nsz, /* Number of nodepoints */
                                          float osz, /* Nodepoints origin */
                                          float dsz /* Nodepoints sampling */)
/*< Use NIP sources location to draw interfaces 
Note: If the velocity model is correct the NIP sources location coincides with interfaces. So, they can be used to draw interface through cubic spline interpolation.
>*/
{
        int nxs=ns; // Number of NIP sources for each interface
        int nxsz=nsz; // Number of nodepoints for each interface
        int i, im; // Loop counter
        float *tsx, *tsz; // Temporary spline vector
        float *coef; // Coefficients matrix
        float xx, xs; // X coordinate
        float oxs; // Spline's origin
        int l; // Spline index

        tsx = sf_floatalloc(nxs+2);
        tsz = sf_floatalloc(nxs+2);
        coef = sf_floatalloc(4*(nxs+2-1));
        for(i=0;i<1;i++){
                for(im=0;im<nxs;im++){
                        tsz[im]=s[im][0];
                        tsx[im]=s[im][1];
			//sf_warning("tsz=%f",tsz[im]);
                }
                sortingXinAscendingOrder(tsx,tsz,nxs);
                calculateSplineCoeficients(nxs,tsx,tsz,coef);
                oxs=tsx[0];
                for(im=0;im<nxsz;im++){
                        xx=im*dsz+osz;
                        if(xx<tsx[0]){
				sz[im]=tsz[0];
			}else if(xx>tsx[nxs-1]){
				sz[im]=tsz[nxs-1];
			}else{
                                l = binarySearch(xx,tsx,nxs);
                                oxs=tsx[l];
                                xs=xx-oxs;
                                sz[im]=coef[l*4+0]*xs*xs*xs+coef[l*4+1]*xs*xs+coef[l*4+2]*xs+coef[l*4+3];
                        }
                }
        }

	//dumpfloat1("tx",tsx,nxs);
	//dumpfloat1("tz",tsz,nxs);
	//dumpfloat1("sz",sz,nxsz);
	free(tsx);
	free(tsz);
	free(coef);
}


void updateVelocityModel3(
			   float *vel,
			   int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			   float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
			   float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
			   float *sv, /* Velocity model disturbance */
			   float *sz,
			   int nsz,
			   float osz,
			   float dsz,
			   bool first)
/*< Velocity model update
Note: This function uses a sv (layers velocity) vector and sz (depth interfaces
coordinates) vector to build the depth velocity model. There is nsv constant
velocity layers in the model and nsv-1 interfaces separating them.
These interfaces are described with nsz control points (nodes) in the sz vector and
they are interpolated using natural cubic spline interpolation.
 >*/
{
        int i, j, i1, i2, l=0;
	float *x=NULL, **coef=NULL;
	float zi[1];
	float xx, z;

	if(first){
		for(i2=0;i2<n[1];i2++){

			for(i1=0;i1<n[0];i1++){
					vel[i2*n[0]+i1]=sv[0];
			}
		}
	}else{

		x = sf_floatalloc(nsz);

		for(i=0;i<nsz;i++)
			x[i] = i*dsz+osz;

	       /* Calculate coefficients matrix (interfaces interpolation) */
		coef = sf_floatalloc2(4*(nsz-1),1);
		calculateSplineCoeficients(nsz,x,sz,coef[0]);

		//zi[0] = (n[0]-1)*d[0]+o[0];

		/* Calculate velocity function */
		for(j=0;j<n[1];j++){

			xx = d[1]*j+o[1];
			if(xx>x[l+1]) l++;
			/* Calculate interfaces z coordinates */
			calcInterfacesZcoord(zi,1,xx-x[l],l,coef);
			for(i=0;i<n[0];i++){
				z = i*d[0]+o[0];
				if(z>zi[0]){
					vel[n[0]*j+i] = sv[0];
				}else{
					vel[n[0]*j+i] = sqrt(1./vel[n[0]*j+i]);
				}
			} /* Loop over depth */
		} /* Loop over distance */

		//dumpfloat1("v",vel,n[0]);
		free(coef);
		free(x);
	}
}

void buildSlownessModelFromVelocityModel3(
					 float *vel,
					 int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
			 		 float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
					 float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
					 float *sv, /* Velociy disturbance */
					 float *sz,
					 int nsz,
					 float osz,
					 float dsz,
					 bool first)
/*< Slowness model build from velocity model
Note: This function is a function wrapper to updateVelocityModel function.
It calls that function to update the velocity model and build the slowness
model matrix using the slowness definition slow=(1.0/(v*v)). 
 >*/
{

	int i, nm; // Loop counters and indexes

	nm =n[0]*n[1];
	updateVelocityModel3(vel,n,o,d,sv,sz,nsz,osz,dsz,first);

	/* transform velocity to slowness */
	for(i=0;i<nm;i++){
			vel[i] = 1.0/(vel[i]*vel[i]);
	}
}

