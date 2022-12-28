/*
* test_raytrace.c (C)
* 
* Purpose: Unit tests of raytrace.c lib.
* 
* Site: https://dirack.github.io
* 
* Version 1.0
* 
* Programmer: Rodolfo A C Neves 26/12/2021
* 
* Email: rodolfo_profissional@hotmail.com
* 
* License: GPL-3.0 <https://www.gnu.org/licenses/gpl-3.0.txt>.
*/

#include "Unity/unity.h"
#include "raytrace.h"
#include <stdio.h>
#include <rsf.h>

void setUp(){};

void tearDown(){};

void test_getTransmissionAngleSnellsLaw()
/*< Function to get transmission angle using snell's law >*/
{
	float ei=SF_PI/6., vt=1.5, vi=1.7;
	TEST_ASSERT_FLOAT_WITHIN(0.001,0.456,getTransmissionAngleSnellsLaw(ei,vi,vt));
}

void test_getTransmitedRNIPHubralTransmissionLaw()
/*< Test function to get Transmited RNIP parameter through interface using Hubral's transmission law >*/
{
	float rnip=1.0;
	float vt=1.5, vi=1.7;
	float ei=SF_PI/6.;
	itf2d it2;
	int i;
	float *sz;
	float kf;

	sz = sf_floatalloc(5);

	for(i=0;i<5;i++)
		sz[i]=1.0;

	it2 = itf2d_init(sz,5,0.,0.5);

	kf = calculateInterfaceCurvature(it2,1.2);
	getTransmitedRNIPHubralTransmissionLaw(&rnip,vt,vi,ei,kf);
	TEST_ASSERT_FLOAT_WITHIN(0.01,1.218,rnip);
}

void test_calculateIncidentAngle()
/*< Test function to get incident angle for a given ray sample coordinate in ray trajectory >*/
{
	float t[4]={1.22,1.05,1.2,1.0};
	float **traj;
	itf2d it2;
	float *sz;
	int i;

	sz = sf_floatalloc(5);
	for(i=0;i<5;i++)
		sz[i]=1.0;

	it2 = itf2d_init(sz,5,0.,0.5);

	traj = sf_floatalloc2(2,2);

	traj[0][0]=t[0];
	traj[0][1]=t[1];
	traj[1][0]=t[2];
	traj[1][1]=t[3];

	TEST_ASSERT_FLOAT_WITHIN(0.01,1.19,calculateIncidentAngle(it2,traj,1));
}

void test_getVelocityForRaySampleLocation()
/*< Test of the function to get velocity from grid for a ray sample location given >*/
{
	float *slow;
	int im;
	raytrace rt;
	int n[2]={10,10};
	float o[2]={0.,-2.};
	float d[2]={0.01,0.01};
	float t[2]={0.05,-1.95};
	float **traj;

        slow =  sf_floatalloc(100);
	traj = sf_floatalloc2(2,1);
	traj[0][0] = t[0];
	traj[0][1] = t[1];

	for(im=0;im<100;im++)
		slow[im] = 1./(1.5*1.5);

	rt = raytrace_init(2,true,10,1.0,n,o,d,slow,4);

	TEST_ASSERT_FLOAT_WITHIN(0.01,1.5,getVelocityForRaySampleLocation(rt,traj,0));
}

float cubicInterface(float x)
/*< Cubic function for tests >*/
{
	return x*x*x-3*x*x+4;
}

void test_firstDerivativeFunction()
/*< Test first derivative numerical calculation >*/
{
	float fx[5];
	float dfdx[5];
	float x;
	float ox=-3.14159;
	float dx=0.01;
	int nx=600;
	int ix;

	for(ix=0;ix<nx;ix++){
		x = ix*dx+ox;
		fx[0]=sinf(x-2*0.01);
		fx[1]=sinf(x-1*0.01);
		fx[2]=sinf(x);
		fx[3]=sinf(x+1*0.01);
		fx[4]=sinf(x+2*0.01);
		first_deriv(0.01,fx,dfdx,5);
		TEST_ASSERT_FLOAT_WITHIN(0.01,cosf(x),dfdx[2]);
	}
}

void test_secondDerivativeFunction()
/*< Test second derivative numerical calculation >*/
{

	TEST_IGNORE_MESSAGE("Second derivative is not exact!");
	float fx[5];
	float dfdx[5];
	float x;
	float ox=-3.14159;
	float dx=0.01;
	int nx=600;
	int ix;
	int j;

	for(ix=0;ix<nx;ix++){
		x = ix*dx+ox;
		for(j=0;j<5;j++)
			fx[j]=sinf(x+(j-2)*0.01);
		second_deriv(0.01,fx,dfdx,5);
		TEST_ASSERT_FLOAT_WITHIN(0.01,-sinf(x),dfdx[2]);
	}
}


void test_calculateInterfaceCurvature()
/*< TODO >*/
{
	TEST_IGNORE_MESSAGE("Curvature calcutation depends on second derivative");
	float capa[4]={6.,0.,6.,0.016};
	itf2d it2;
	int i;
	float x[5]={-1,0,1,2,3};
	float y[10]={0.,3.125,4.0,3.375,2.,0.625,0.,0.875,4.,10.125};

	it2 = itf2d_init(y,10,-1,0.5);

	for(i=0;i<4;i++)
		TEST_ASSERT_FLOAT_WITHIN(0.1,capa[i],calculateInterfaceCurvature(it2,i));
}

void test_sortingXinAscendingOrder()
/*TODO*/
{
	int i;
	float x[13]={4.27705002, 3.03945422, 3.26645589, 2.03656244, 2.27515054, 2.51043391, 2.78882909, 4.02767229, 1.50065649, 1.79396605,1.79396605, 3.51551676, 3.76459265};
	float y[13]={0.966627121, 0.966370642, 0.963363767, 0.966482043, 0.974164248, 0.969641984, 0.972658157, 0.959083617, 0.972544372, 0.974167109, 0.974167109, 0.963503361, 0.969640553};
	float xs[13]={1.50065649, 1.79396605, 1.79396605, 2.03656244, 2.27515054, 2.51043391, 2.78882909, 3.03945422, 3.26645589, 3.51551676,3.76459265, 4.02767229, 4.27705002};
	float ys[13]={0.972544372, 0.974167109, 0.974167109, 0.966482043, 0.974164248, 0.969641984, 0.972658157, 0.966370642, 0.963363767,0.963503361, 0.969640553, 0.959083617, 0.966627121};
	sortingXinAscendingOrder(x,y,13);
	for(i=1;i<13;i++){
		TEST_ASSERT_FLOAT_WITHIN(0.01,xs[i],x[i]);
		TEST_ASSERT_FLOAT_WITHIN(0.01,ys[i],y[i]);
	}
}

void test_binarySearch()
/**/
{
	float x[10]={0.,1.,2.,3.,4.,5.,6.,7.,8.,9.};
	TEST_ASSERT_EQUAL(4,binarySearch(4.5,x,10));
	TEST_ASSERT_EQUAL(4,binarySearch(4.001,x,10));
	TEST_ASSERT_EQUAL(4,binarySearch(4.999,x,10));
	TEST_ASSERT_EQUAL(6,binarySearch(6.5,x,10));
	TEST_ASSERT_EQUAL(6,binarySearch(6.001,x,10));
	TEST_ASSERT_EQUAL(6,binarySearch(6.999,x,10));
	TEST_ASSERT_EQUAL(0,binarySearch(0.5,x,10));
	TEST_ASSERT_EQUAL(0,binarySearch(0.001,x,10));
	TEST_ASSERT_EQUAL(0,binarySearch(0.999,x,10));
	TEST_ASSERT_EQUAL(8,binarySearch(8.5,x,10));
	TEST_ASSERT_EQUAL(8,binarySearch(8.001,x,10));
	TEST_ASSERT_EQUAL(8,binarySearch(8.999,x,10));
}

int main(void){

	UNITY_BEGIN();
	/*RUN_TEST(test_getTransmissionAngleSnellsLaw);
	RUN_TEST(test_getTransmitedRNIPHubralTransmissionLaw);
	RUN_TEST(test_calculateIncidentAngle);
	RUN_TEST(test_getVelocityForRaySampleLocation);*/
	RUN_TEST(test_firstDerivativeFunction);
	RUN_TEST(test_secondDerivativeFunction);
	RUN_TEST(test_calculateInterfaceCurvature);
	RUN_TEST(test_sortingXinAscendingOrder);
	RUN_TEST(test_binarySearch);
	return UNITY_END();
}
