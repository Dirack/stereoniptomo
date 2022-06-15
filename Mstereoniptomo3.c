/* Landa 1988 experiment: VFSA velocity inversion based on stereotomography and NIP tomography strategies

This program is a reproduction of the experiment in the article 'A method for determination of velocity and depth from seismic reflection data' from Landa 1988.

That program do the model setup using parameters (t0,m0,BETA,RNIP) picked from stacked section. Rays are traced into the model from m0 and using BETA as the initial direction until t0/2 traveltime is consumed. The final ray coordinate is the NIP source position.

So, the forward modeling is done by ray tracying, from NIP sources to acquisition surface. The program gets reflection traveltime for a set of reflection ray pairs.

The semblance is calculated stacking prestack data amplitudes using Non-Hyperbolic CRS traveltime approximation in RN=RNIP limit (CDS condition). This semblance is used as a convergence criteria for VFSA global optimization algorithm to obtain optimized velocity model. This is done layer by layer.

*/

#define DEBUG_UTILS
#ifdef DEBUG_UTILS
#include "utils.h"
#endif
#include <math.h>
#include <rsf.h>
#include <time.h>
#include "tomography.h"
#include "vfsacrsnh_lib.h"
#include "velocity_lib.h"

#define RAD2DEG 180./SF_PI

int main(int argc, char* argv[])
{
	bool verb; // Verbose parameter
	int n[2]; // Velocity grid dimensions n[0]=n1, n[1]=n2
	float d[2]; // Velocity grid sampling d[0]=d1, d[1]=d2
	float o[2]; // Velocity grid origin o[0]=o1, o[1]=o2
	float** s; // NIP sources position (z,x)
	float cnewv[1]; // Temporary parameters vector used in VFSA
	float otsv[1]; // Optimized parameters vector
	float tmis0; // Best time misfit
	float otmis=0; // Best time misfit
	float deltaE; // Delta (Metrópolis criteria in VFSA)
	float Em0=0; // Energy (VFSA algorithm)
	float PM; // Metrópolis criteria
	float temp=1; // Temperature for VFSA algorithm
	float u=0; // Random number between 0 and 1
	int nit; // Number of VFSA iterations
	float temp0; // Initial temperature for VFSA
	float c0; // Damping factor for VFSA
	int ndim; // n1 dimension in shotsfile, should be equal 2
	int nshot; // n2 dimensions in shotsfile, number of shots
	int nm; // Number of samples in velocity grid n1*n2
	float* a; // Normal Ray initial angle for each NIP source
	float* slow; // slowness model
	int im; // loop counter
	float v; // Velocity temporary variable
	float v0; // Near surface velocity
	int ns; // Number of NIP sources
	int q; // Loop counter for VFSA iteration
	float tmis; // data time misfit value
	float *m0; // CMP's for normal rays
	float *t0; // t0's for normal rays
	float *RNIP; // Rnip parameters vector
	float *BETA; // Beta parameters vector
	float sv[1]; // Layer's Velocity
	float minvel; // Minimun layer velocity
	float maxvel; // Maximum layer velocity
	bool first; // First interface inversion?
	bool base; // Last interface inversion?
	float **ots; // NIP sources from inversion
	float *otsz; // Interface Z spline nodes (optimal)
	float *sz; // Interface Z spline nodes
	int nsz; // Number of interfaces nodes
	float osz; // Interfaces nodes origin
	float dsz; // Interfaces nodes sampling
	float ***data; // Prestack data A(m,h,t)
	int data_n[3]; // n1, n2, n3 dimension of data
	float data_o[3]; // o1, o2, o3 axis origins of data
	float data_d[3]; // d1, d2, d3 sampling of data
	sf_file shots; // NIP sources (z,x)
	sf_file vel; // background velocity model
	sf_file vz_file; // Initial Layer velocity
	sf_file sz_file; // Initial interface splines nodes
	sf_file zspline; // Final interface splines nodes
	sf_file velinv; // Inverted velocity model
	sf_file m0s; // Central CMPs m0
	sf_file t0s; // Normal ray traveltimes
	sf_file rnips; // RNIP parameter for each m0
	sf_file betas; // BETA parameter for each m0
	sf_file vspline; // Layers velocity (output)
	sf_file datafile; // Prestack data A(m,h,t)

	sf_init(argc,argv);

	shots = sf_input("shotsfile");
	vel = sf_input("in");
	vz_file = sf_input("sv");
	sz_file = sf_input("sz");
	velinv = sf_output("out");
	vspline = sf_output("vspline");
	zspline = sf_output("zspline");
	m0s = sf_input("m0s");
	t0s = sf_input("t0s");
	rnips = sf_input("rnips");
	betas = sf_input("betas");
	datafile = sf_input("data");

	/* Velocity model: get 2D grid parameters */
	if(!sf_histint(vel,"n1",n)) sf_error("No n1= in input");
	if(!sf_histint(vel,"n2",n+1)) sf_error("No n2= in input");
	if(!sf_histfloat(vel,"d1",d)) sf_error("No d1= in input");
	if(!sf_histfloat(vel,"d2",d+1)) sf_error("No d2= in input");
	if(!sf_histfloat(vel,"o1",o)) o[0]=0.;
	if(!sf_histfloat(vel,"o2",o+1)) o[1]=0.;

	/* Prestack data cube */
	if(!sf_histint(datafile,"n1",data_n)) sf_error("No n1= in data");
	if(!sf_histint(datafile,"n2",data_n+1)) sf_error("No n2= in data");
	if(!sf_histint(datafile,"n3",data_n+2)) sf_error("No n3= in data");
	if(!sf_histfloat(datafile,"d1",data_d)) sf_error("No d1= in data");
	if(!sf_histfloat(datafile,"d2",data_d+1)) sf_error("No d2= in data");
	if(!sf_histfloat(datafile,"d3",data_d+2)) sf_error("No d3= in data");
	if(!sf_histfloat(datafile,"o1",data_o)) sf_error("No o1= in data");
	if(!sf_histfloat(datafile,"o2",data_o+1)) sf_error("No o2= in data");
	if(!sf_histfloat(datafile,"o3",data_o+2)) sf_error("No o3= in data");
	
	if(!sf_getbool("verb",&verb)) verb=true;
	/* verbose parameter (y/n) */

	if(!sf_getfloat("v0",&v0)) v0=1.5;
	/* Near surface velocity (Km/s) */

	if(!sf_getint("nit",&nit)) nit=1;
	/* Number of VFSA iterations */

	if(!sf_getfloat("temp0",&temp0)) temp0=5;
	/* Initial temperature for VFSA algorithm */

	if(!sf_getfloat("c0",&c0)) c0=0.1;
	/* Damping factor for VFSA algorithm */

	if(!sf_getfloat("minvel",&minvel)) minvel=1.5;
	/* Layers minimum velocity */

	if(!sf_getfloat("maxvel",&maxvel)) maxvel=2.0;
	/* Layers maximum velocity */

	if(!sf_getbool("first",&first)) first=true;
	/* First interface inversion (y/n) */

	if(!sf_getbool("base",&base)) base=false;
	/* No interface inversion, just output the velocity model (y/n) */

	/* Shotsfile: get shot points (z=0, x=m0) */
	if(!sf_histint(shots,"n1",&ndim) || 2 != ndim)
		sf_error("Must have n1=2 in shotsfile");
	if(!sf_histint(shots,"n2",&nshot)) sf_error("No n2= in shotsfile");
	s = sf_floatalloc2(ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	sf_fileclose(shots);
	ns=nshot;

	if(!sf_histint(sz_file,"n1",&nsz)) sf_error("No n1= in sz_file");
	if(!sf_histfloat(sz_file,"o1",&osz)) sf_error("No o1= in sz_file");
	if(!sf_histfloat(sz_file,"d1",&dsz)) sf_error("No d1= in sz_file");
	sz = sf_floatalloc(nsz);
	sf_floatread(sz,nsz,sz_file);
	otsz = sf_floatalloc(nsz);
	ots = sf_floatalloc2(ndim,ns);

	/* Read initial velocity from file */
	sf_floatread(sv,1,vz_file);

	/* Read prestack data cube A(t,h,m) */
	data = sf_floatalloc3(data_n[0],data_n[1],data_n[2]);
	sf_floatread(data[0][0],data_n[0]*data_n[1]*data_n[2],datafile);

	/* allocate parameters vectors */
	a = sf_floatalloc(ns);
	m0 = sf_floatalloc(ns);
	sf_floatread(m0,ns,m0s);
	t0 = sf_floatalloc(ns);
	sf_floatread(t0,ns,t0s);
	RNIP = sf_floatalloc(ns);
	sf_floatread(RNIP,ns,rnips);
	BETA = sf_floatalloc(ns);
	sf_floatread(BETA,ns,betas);

	/* get slowness squared (Background model) */
	nm = n[0]*n[1];
	slow =  sf_floatalloc(nm);
	sf_floatread(slow,nm,vel);

	for(im=0;im<nm;im++){
		v = slow[im];
		slow[im] = 1./(v*v);
	}

	if(verb){
		sf_warning("Command line Parameters");
		sf_warning("v0=%f nit=%d temp0=%f c0=%f",v0,nit,temp0,c0);
		sf_warning("Input file (Velocity model)");
		sf_warning("n1=%d d1=%f o1=%f",*n,*d,*o);
		sf_warning("n2=%d d2=%f o2=%f",*(n+1),*(d+1),*(o+1));
		sf_warning("Input file (Prestack data)");
		sf_warning("n1=%d d1=%f o1=%f",*data_n,*data_d,*data_o);
		sf_warning("n2=%d d2=%f o2=%f",*(data_n+1),*(data_d+1),*(data_o+1));
		sf_warning("n3=%d d3=%f o3=%f",*(data_n+2),*(data_d+2),*(data_o+2));
		sf_warning("Input file (shotsfile)");
		sf_warning("n1=%d",ndim);
		sf_warning("n2=%d",nshot);
		sf_warning("Input file (anglefile, t0s, m0s, rnips, betas)");
		sf_warning("n1=%d",ns);
		sf_warning("Input file (vz) - Initial velocity");
		sf_warning("vi=%f",sv[0]);
	}

	/* Velocity model from inversion */
	sf_putint(velinv,"n1",n[0]);
	sf_putint(velinv,"n2",n[1]);
	sf_putint(velinv,"n3",1);
	sf_putfloat(velinv,"d1",d[0]);
	sf_putfloat(velinv,"d2",d[1]);
	sf_putfloat(velinv,"o1",o[0]);
	sf_putfloat(velinv,"o2",o[1]);
	sf_putfloat(velinv,"d3",1);
	sf_putfloat(velinv,"o3",0);

	/* velocity and interfaces (output) */
	sf_putint(vspline,"n1",1);
	sf_putint(vspline,"n2",1);

	/* Build initial velocity model and setup NIP sources */
	buildSlownessModelFromVelocityModel3(slow,n,o,d,sv,sz,nsz,osz,dsz,first);
	modelSetup(s, ns,  m0, t0, BETA,  a,  n,  d,  o,  slow);
	tmis0=0.;
	otmis=tmis0;

	/* Initiate optimal parameters vectors */
	for(im=0;im<ns;im++){
		ots[im][0]=s[im][0];
		ots[im][1]=s[im][1];
	}
	otsv[0]=sv[0];

	srand(time(NULL));


	/* Very Fast Simulated Annealing (VFSA) algorithm */
	if(!base){
		for (q=0; q<nit; q++){
		
			/* calculate VFSA temperature for this iteration */
			temp=getVfsaIterationTemperature(q,c0,temp0);
							
			/* parameter disturbance */
			disturbParameters3(temp,cnewv,sv,minvel,maxvel,1);

			/* Update velocity model */
			buildSlownessModelFromVelocityModel3(slow,n,o,d,cnewv,sz,nsz,osz,dsz,first);

			/* NIP sources setup for new model */
			tmis=0;
			modelSetup(s, ns,  m0, t0, BETA,  a,  n,  d,  o,  slow);

			/* Forward modeling */
			// TODO change tmis variable name to semb (Semblance)
			tmis=forwardModeling(s,v0,t0,m0,RNIP,BETA,n,o,d,slow,a,ns,data,data_n,data_o,data_d);
		
			if(fabs(tmis) > fabs(tmis0) ){
				otmis = fabs(tmis);
				/* Keep optimized parameters */
				for(im=0;im<ns;im++){
					ots[im][0]=s[im][0];
					ots[im][1]=s[im][1];
				}
				otsv[0]=cnewv[0];
				tmis0 = fabs(tmis);
			}

			/* VFSA parameters update condition */
			deltaE = fabs(tmis) - Em0;
			
			/* Metrópolis criteria */
			PM = expf(-deltaE/temp);
			
			if (deltaE<=0){
				sv[0]=cnewv[0];
				Em0 = fabs(tmis);
			} else {
				u=getRandomNumberBetween0and1();
				if (PM > u){
					sv[0]=cnewv[0];
					Em0 = fabs(tmis);
				}
			}	
			
			// TODO do not show this iteration semblance, only optimal one	
			sf_warning("%d/%d Missfit(%f) vel=%f v=%f %f ;",q+1,nit,otmis,otsv[0],cnewv[0],tmis);

		} /* loop over VFSA iterations */
	} /* If base=true skip VFSA, only generate the velocity model */

	/* Generate optimal velocity model and interfaces */
	updateVelocityModel3(slow,n,o,d,otsv,sz,nsz,osz,dsz,first);

	interfaceInterpolationFromNipSources(ots,ns,otsz,nsz,osz,dsz);

	/* Write interfaces nodepoints */
	sf_putint(zspline,"n1",nsz);
	sf_putint(zspline,"n2",1);
	sf_putfloat(zspline,"o1",osz);
	sf_putfloat(zspline,"d1",dsz);

	/* Write layer velocity */
	sf_floatwrite(otsv,1,vspline);
	sf_floatwrite(otsz,nsz,zspline);
	sf_floatwrite(slow,nm,velinv);
}
