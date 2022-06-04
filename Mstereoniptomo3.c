/* Landa 1988 experiment: VFSA velocity inversion based on stereotomography and NIP tomography strategies

This program is a reproduction of the experiment in the article 'A method for determination of velocity and depth from seismic reflection data' avaliable in the doc directory of this repository

The initial velocity model and NIP sources position used in this program is set up using sfnipmodsetup, that program does the forward modeling by ray tracying, from NIP sources to acquisition surface and gets reflection traveltime for a set of reflection ray pairs.

The time misfit is calculated by the difference between the reflection traveltime obtained in the forward modeling and the traveltime calculated by CRE traveltime approximation formula, for RNIP and BETA parameters given. This time misfit is used as a convergence criteria for VFSA global optimization algorithm to obtain optimized velocity model.

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
	float minvel;
	float maxvel;
	float ***data; // Prestack data A(m,h,t)
	int data_n[3]; // n1, n2, n3 dimension of data
	float data_o[3]; // o1, o2, o3 axis origins of data
	float data_d[3]; // d1, d2, d3 sampling of data
	sf_file shots; // NIP sources (z,x)
	sf_file vel; // background velocity model
	sf_file vz_file;
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
	velinv = sf_output("out");
	vspline = sf_output("vspline");
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

	/* Shotsfile: get shot points */
	if(!sf_histint(shots,"n1",&ndim) || 2 != ndim)
		sf_error("Must have n1=2 in shotsfile");
	if(!sf_histint(shots,"n2",&nshot)) sf_error("No n2= in shotfile");
	s = sf_floatalloc2(ndim,nshot);
	sf_floatread(s[0],ndim*nshot,shots);
	sf_fileclose(shots);
	ns=nshot;

	/* Cubic spline vector */
	/*if(!sf_histint(vz_file,"n1",nsv)) sf_error("No n1= in sv file");
	if(!sf_histint(vz_file,"n2",nsv+1)) sf_error("No n1= in sv file");
	if(!sf_histfloat(vz_file,"o1",osv)) sf_error("No o2= in sv file");
	if(!sf_histfloat(vz_file,"o2",osv+1)) sf_error("No o3= in sv file");
	if(!sf_histfloat(vz_file,"d1",dsv)) sf_error("No d2= in sv file");
	if(!sf_histfloat(vz_file,"d2",dsv+1)) sf_error("No d3= in sv file");*/

	/* Build cubic spline velocity matrix */
	//sv = sf_floatalloc(1);
	sf_floatread(sv,1,vz_file);

	/* VFSA parameters vectors */
	//cnewv = sf_floatalloc(1);
	//otsv = sf_floatalloc(1);

	/* Read prestack data cube */
	data = sf_floatalloc3(data_n[0],data_n[1],data_n[2]);
	sf_floatread(data[0][0],data_n[0]*data_n[1]*data_n[2],datafile);

	a = sf_floatalloc(ns);

	/* allocate parameters vectors */
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
		sf_warning("Input file (vz)");
		//sf_warning("n1=%d d1=%f o1=%f",*nsv,*dsv,*osv);
		//sf_warning("n2=%d d2=%f o2=%f",*(nsv+1),*(dsv+1),*(osv+1));
	}

	/* Use previous misfit as the initial misfit value */
	buildSlownessModelFromVelocityModel3(slow,n,o,d,sv);
	modelSetup(s, ns,  m0, t0, BETA,  a,  n,  d,  o,  slow);
	tmis0=0;//forwardModeling(s,v0,t0,m0,RNIP,BETA,n,o,d,slow,a,ns,data,data_n,data_o,data_d);
	otmis=tmis0;

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
	/*sf_putfloat(vspline,"o1",osv[0]);
	sf_putfloat(vspline,"o2",osv[1]);
	sf_putfloat(vspline,"d1",dsv[0]);
	sf_putfloat(vspline,"d2",dsv[1]);*/
	
	/* Intiate optimal parameters vectors */
	//for(im=0;im<nsv[0]*nsv[1];im++)
	otsv[0]=sv[0];

	srand(time(NULL));
	/* Very Fast Simulated Annealing (VFSA) algorithm */
	for (q=0; q<nit; q++){
	
		/* calculate VFSA temperature for this iteration */
		temp=getVfsaIterationTemperature(q,c0,temp0);
						
		/* parameter disturbance */
		disturbParameters3(temp,cnewv,sv,minvel,maxvel,1);

	//	dumpfloat1("cnewv",cnewv,nsv[0]*nsv[1]);

		/* Function to update velocity model */
		buildSlownessModelFromVelocityModel3(slow,n,o,d,cnewv);

		//sf_warning("%d",__LINE__);
		tmis=0;
		modelSetup(s, ns,  m0, t0, BETA,  a,  n,  d,  o,  slow);
		//sf_warning("%d",__LINE__);
		//dumpfloat1("slow",slow,n[[0]]);
		tmis=forwardModeling(s,v0,t0,m0,RNIP,BETA,n,o,d,slow,a,ns,data,data_n,data_o,data_d);
		//dumpfloat1("slow",slow,nm);
		//sf_error("capa");

	
		if(fabs(tmis) > fabs(tmis0) ){
			otmis = fabs(tmis);
			/* optimized parameters */
			//for(im=0;im<nsv[0]*nsv[1];im++)
			otsv[0]=cnewv[0];
			tmis0 = fabs(tmis);
		}

		/* VFSA parameters update condition */
		deltaE = fabs(tmis) - Em0;
		
		/* Metrópolis criteria */
		PM = expf(-deltaE/temp);
		
		if (deltaE<=0){
			//for(im=0;im<nsv[0]*nsv[1];im++)
			sv[0]=cnewv[0];
			Em0 = fabs(tmis);
		} else {
			u=getRandomNumberBetween0and1();
			if (PM > u){
				//for(im=0;im<nsv[0]*nsv[1];im++)
				sv[0]=cnewv[0];
				Em0 = fabs(tmis);
			}
		}	
			
		sf_warning("%d/%d Missfit(%f) ;",q+1,nit,otmis);

	} /* loop over VFSA iterations */

	/* Generate optimal velocity model */
	updateVelocityModel3(slow,n,o,d,otsv);

	/* Write velocity cubic spline function */
	sf_floatwrite(otsv,1,vspline);
	sf_floatwrite(slow,nm,velinv);
}
