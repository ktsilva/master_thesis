/*
 *  GSLmyLib.h - to help me work with gsl faster!!!! vrum vrum
 *  all start with "cat"
 *
 *  Created by Catarina Silva Fernandes 2010.
 *
 */

#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"
//large collection of random number generators which can be accessed through a uniform interface

//compile w/ gsl: c++ -Wall -I/local/opt/star/include -L/local/opt/star/lib -l gsl


// *******  INIT RANDOM SEED ET AL.  ******* //
void cat_rng_init(){
	srand(time(NULL)); //generates seed æhæhæh 
	
}

// *******  UNIFORM DISTRIBUTION  ******* //

//returns a rdm number (double) uniformily between [0,1) - default by GNU
double cat_rng_uniform_double_zero_one(){
	const gsl_rng_type * T; //holds static information about each type of generator
	gsl_rng * r; //describes an instance of a generator created from a given gsl_rng_type
	gsl_rng_env_setup(); /*reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED 
						  and uses their values to set the corresponding library variables 
						  gsl_rng_default and gsl_rng_default_seed.*/
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T); //returns a pointer to a newly-created instance of a random number generator of type T
	gsl_rng_set(r, rand()); //initializes (or `seeds') the random number generator. it will be always different because we are generating a seed every time :)
	
	double uni;
	uni = gsl_rng_uniform (r);
	
	return uni;

}

//returns a rdm number (double) uniformily between [min,max) 
double cat_rng_uniform_double_min_max(double min, double max) {
	
	const gsl_rng_type * T; //holds static information about each type of generator
	gsl_rng * r; //describes an instance of a generator created from a given gsl_rng_type
	gsl_rng_env_setup(); /*reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED 
						  and uses their values to set the corresponding library variables 
						  gsl_rng_default and gsl_rng_default_seed.*/
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T); //returns a pointer to a newly-created instance of a random number generator of type T
	gsl_rng_set(r, rand()); //initializes (or `seeds') the random number generator. it will be always different because we are generating a seed every time :)
	
	double uni;
	uni = min + (max-min)*gsl_rng_uniform (r);
	
	return uni;

}

//returns a rdm number (int) uniformily between [0,max) - default by GNU
int cat_rng_uniform_int_zero_max(int max){
	const gsl_rng_type * T; //holds static information about each type of generator
	gsl_rng * r; //describes an instance of a generator created from a given gsl_rng_type
	gsl_rng_env_setup(); /*reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED 
						  and uses their values to set the corresponding library variables 
						  gsl_rng_default and gsl_rng_default_seed.*/
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T); //returns a pointer to a newly-created instance of a random number generator of type T
	gsl_rng_set(r, rand()); //initializes (or `seeds') the random number generator. it will be always different because we are generating a seed every time :)
	
	int uni;
	uni = gsl_rng_uniform_int (r, max);
		
	//If n is larger than the range of the generator then the function calls the error handler with an error code of GSL_EINVAL and returns zero. 
	return uni;
}

//returns a rdm number (int) uniformily between [min,max)
int cat_rng_uniform_int_min_max(int min, int max){
	const gsl_rng_type * T; //holds static information about each type of generator
	gsl_rng * r; //describes an instance of a generator created from a given gsl_rng_type
	gsl_rng_env_setup(); /*reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED 
						  and uses their values to set the corresponding library variables 
						  gsl_rng_default and gsl_rng_default_seed.*/
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T); //returns a pointer to a newly-created instance of a random number generator of type T
	gsl_rng_set(r, rand()); //initializes (or `seeds') the random number generator. it will be always different because we are generating a seed every time :)
	
	int uni;
	uni = min+ (max-min)*gsl_rng_uniform_int (r, max);
	return uni;
}


//returns a rdm number (double) under a gaussian disribution
double cat_rng_gaussian(double sigma) {
	
	const gsl_rng_type * T; //holds static information about each type of generator
	gsl_rng * r; //describes an instance of a generator created from a given gsl_rng_type
	gsl_rng_env_setup(); /*reads the environment variables GSL_RNG_TYPE and GSL_RNG_SEED 
						  and uses their values to set the corresponding library variables 
						  gsl_rng_default and gsl_rng_default_seed.*/
	
	T = gsl_rng_default;
	r = gsl_rng_alloc (T); //returns a pointer to a newly-created instance of a random number generator of type T
	gsl_rng_set(r, rand()); //initializes (or `seeds') the random number generator. it will be always different because we are generating a seed every time :)
	
	double uni;
	uni = gsl_ran_gaussian (r, sigma);
	
	return uni;
	
}

