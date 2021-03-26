/*
 *  MyProgramJumpMteCarlo.h
 *  
 *
 *  Created by Catarina Silva Fernandes on2010.
 * 
 *
 */

#include <stdlib.h>
#include <time.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <locale>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include "MyProgramTestingFit.h"
#include "libraries/arrayErrorPropagation.h"

#include <time.h>


array * jeans_mass_from_total_mass(array TotalMass, array GasMass, array Radius, 
								   array Disp_velocity, array error_Disp_velocity,
								   array beta_profile, array error_beta_profile, double * guess) {
	
	int nBin = TotalMass.numberBins();
	
	array * jeans_mass; jeans_mass = new array[2];
	jeans_mass[0].init(nBin); jeans_mass[1].init(nBin);
	
	array dm_density, error_dm_density, dm_mass, error_dm_mass;
	dm_density.init(nBin); dm_mass.init(nBin);
	error_dm_density.init(nBin); error_dm_mass.init(nBin);
	
	dm_mass = ArrayMinusArray(TotalMass, GasMass);
	
	array * fitting; fitting = new array[2];
	fitting = fit_mass_from_density_PL(Radius, dm_mass, guess); //return fit for dm_mass & (parameter +/- ∆ parameter)
	printf("DM density profile - rho0 = %g +/- %g; r_s = %g +/- %g; alpha = %g +/- %g; beta = %g +/- %g\n", 
		   fitting[1].read(0), fitting[1].read(4), fitting[1].read(1), fitting[1].read(5),
		   fitting[1].read(2), fitting[1].read(6), fitting[1].read(3), fitting[1].read(7));
	
	//we don't need this
	//error_dm_mass = err_mass_assume_density_powerLaw(fitting, Radius) ;
	//printf("Fitting dm_mass. i.e., m[n/2] = %g +/- %g\n", dm_mass[0.5*nBin], error_dm_mass[0.5*nBin]);
	
	//dm_density = DensityFromMassSpherical(Radius, dm_mass);
	//error_dm_density = error_one_variable(dm_density, dm_mass, error_dm_mass);
	dm_density = DensityPowerLaw(Radius, fitting[1].read(0), fitting[1].read(1), fitting[1].read(2), fitting[1].read(3));
	
	double * fit; fit = new double[8];
	for (int i = 0; i < nBin; i++) fit[i] = fitting[1].read(i);
	
	error_dm_density = err_density_assume_powerLaw(fit, Radius);
	printf("Fitting dm_density. i.e., m[n/2] = %g +/- %g\n", dm_density.read(0.5*nBin), error_dm_density.read(0.5*nBin));
	
	jeans_mass[0] = MassFromJeans(Radius, dm_density, ArraySQ(Disp_velocity), beta_profile);
	jeans_mass[1] = err_three_variables(jeans_mass[0], 
										dm_density, error_dm_density, 
										Disp_velocity, error_Disp_velocity, 
										beta_profile, error_beta_profile);
	
	printf("Arrived to new Jeans mass from Non thermal profile - calculating likelihood\n");
	delete[] fit;
	return jeans_mass;
}
array side1_slope_plateau_xLog(array xAxis, double slope1, double x_ratio_slope_to_plateau){
	
	int nBin = xAxis.numberBins();
	array C; C.init(nBin);
	double value;
	
	for (int i = 0; i < nBin; i++) {
		
		//if (xAxis.read(i) < x_chage_slope_to_plateau) value = slope1*log10(xAxis.read(i))+ x_chage_slope_to_plateau;			
		
		//else value = 0.;
		
		//C.set(i,value);
	}
	return C;
}
//tenho de dar o valor de R200 se nao esta merda nao funciona!!!
//init_conditions[2] = .5*(.8-init_conditions[0]*log10(.01*.001));	//position x0_1
//	init_conditions[3] = .5*(-.8-init_conditions[1]*log10(10.*1.));//* 1.e3*1.8e3*onePcMtr);	//position x0_2

array side2_slope_plateau_xLog(array xAxis, double slope2, double x_chage_plateau_to_slope){
	
	int nBin = xAxis.numberBins();
	array C; C.init(nBin);
	double value;
	
	for (int i = 0; i < nBin; i++) {
		if (xAxis.read(i) < x_chage_plateau_to_slope) value = 0.;			
		
		else value = slope2*log10(xAxis.read(i)) + x_chage_plateau_to_slope;
		
		C.set(i,value);
	}
	return C;
}

double likelihood (array A, array error_A, array B, array error_B){
	
	int nBin = A.numberBins();
	double value = 0.;
	for (int i = 0; i < nBin; i++) {
		//value += pow(A.read(i) - B.read(i),2.)/ (A.read(i) + B.read(i)); 
		value += pow(A.read(i) - B.read(i),2.)/ (pow(error_A.read(i),2.)+ pow(error_B.read(i),2.));
	}
	
	return value/(nBin-4.);
	
}

void metropolis_algorithm_same_delta(string name_file, int nIter, int nRej, double delta, double * x_i, double * min_max_plateau,
									 array rkpc, array radius, array gas_density, array gas_temperature, 
									 array disp_velocity, array dm_beta, double * dm_dens_guess) {
	
	ofstream outFile;
	string MyFile;
	int nBin = rkpc.numberBins();
	
	double r0l, r1l, r0r, r1r, min_plateau, max_plateau;
	r0l = x_i[0]; r1l = x_i[1]; r0r = x_i[2]; r1r = x_i[3];
	min_plateau = min_max_plateau[0]; max_plateau = min_max_plateau[1];
	
	array gas_mass, mass_HE, mass_true, dm_mass, dm_density, mass_jeans, bias;
	gas_mass.init(nBin); mass_HE.init(nBin); mass_true.init(nBin); dm_mass.init(nBin); dm_density.init(nBin); mass_jeans.init(nBin); bias.init(nBin);
	
	gas_mass = massFromDensityAssumeSpherical(radius, gas_density);
	mass_HE = MassFromHE(radius,gas_density,gas_temperature);
	bias = funct_elephant_inside_snake(rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
	//	writeArrayToFile(rkpc, DoubleTimesArray(100.,bias), name_file+"_bias_init");
	
	mass_true = total_mass_from_bias_elephant_inside_snake(mass_HE, rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
	//	writeArrayToFile(rkpc, ArrayDividesDouble(mass_true, oneMsunkg), name_file+"_mass_true_w_ansatz");
	
	dm_mass = ArrayMinusArray(mass_true, gas_mass);
	//fit to DM density
	double * fit_parameters; fit_parameters = new double [8];
	fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	//	writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_w_ansatz");
	
	//	MyFile = name_file + "_chiSQ_r_true_jeans_ansatz.txt";
	//	outFile.open(MyFile.c_str(),ios::out);
	
	double chi_sq_test=0.;
	
	for (int i = 0; i < nBin; i++) {
		double value = pow(mass_true.read(i)-mass_jeans.read(i),2.)/(pow(mass_true.read(i),2.)+ pow(mass_jeans.read(i),2.));
		//		outFile << rkpc.read(i) << " " << value << endl;
		if (rkpc.read(i) >= 0.001 && rkpc.read(i) <= 10.) {
			chi_sq_test+= value;
		}
	}
	//	outFile.close();
	
	
	// INIT METROPOLIS
	
	double point_a,point_b, point_c, point_d, ref;
	double chi_sq[nIter], x0[nIter], x1[nIter], x2[nIter], x3[nIter];
	chi_sq[0] = chi_sq_test;
	x0[0]= r0l; x1[0]= r1l; x2[0]= r0r; x3[0]= r1r;
	
	//	MyFile = name_file + "_chiSQ_vs_X0123.txt";
	//	outFile.open(MyFile.c_str(),ios::out);
	//	outFile << 0 << " " << chi_sq[0] << " " << x0[0] << " " << x1[0] << " " << x2[0] << " " << x3[0]<< endl;
	
	//	MyFile = name_file + "_best_chi_vs_X0123.txt";
	//	ofstream outFile1; outFile1.open(MyFile.c_str(),ios::out);
	
	/*	
	 double tstart, tstop, ttime;
	 
	 tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
	 
	 int i =1; int rej = 0.;
	 while (i<nIter){
	 
	 printf("METROPOLIS: New configuration %d chiSQ = %g \n", i, chi_sq[i-1] );
	 bool make_new = true;
	 while (make_new == true) {
	 
	 r0l = x0[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 r1l = x1[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 r0r = x2[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 r1r = x3[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 
	 // DO we need to set a constraint on r1l and r1r???
	 
	 mass_true = total_mass_from_bias_elephant_inside_snake(mass_HE, rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
	 make_new = false;
	 for (int i = 0; i < nBin; i++) {
	 if (mass_true.read(i) < gas_density.read(i)) make_new = true;
	 }
	 }
	 
	 dm_mass = ArrayMinusArray(mass_true, gas_density);
	 
	 fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	 
	 dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1],
	 fit_parameters[2],fit_parameters[3]);
	 
	 mass_jeans= MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	 
	 chi_sq_test=0.;
	 for (int j = 0; j < nBin; j++) {
	 if (rkpc.read(j) >= 0.001 && rkpc.read(j) <= 10.) {
	 chi_sq_test+= pow(mass_true.read(j)-mass_jeans.read(j),2.)/(pow(mass_true.read(j),2.)+ pow(mass_jeans.read(j),2.));
	 }
	 }
	 //printf("testing chiSq  = %g \n", chi_sq_test);
	 
	 if(chi_sq_test < chi_sq[i-1]){ // accepts
	 chi_sq[i] = chi_sq_test;
	 x0[i]= r0l; x1[i]= r1l; 
	 x2[i]= r0r; x3[i]= r1r;
	 point_a = r0l; point_b = r1l;
	 point_c = r0r; point_d = r1r;
	 if (chi_sq_test < ref){
	 ref = chi_sq_test;
	 outFile1 << i << " " << ref << " " << point_a << " " << point_b << " " << point_c << " " << point_d << endl;
	 }
	 rej = 0.;
	 
	 }
	 else { //there is a small probability of accepting
	 double r = cat_rng_uniform_double_zero_one();
	 double w = exp(-.5*(chi_sq_test - chi_sq[i-1]));
	 //	printf("i=%d r=%g  w=%g\n", i, r, w );
	 if (r > w) { // i do not accept the new configuration
	 chi_sq[i] = chi_sq[i-1];
	 x0[i]= x0[i-1]; x1[i]= x1[i-1]; x2[i]= x2[i-1]; x3[i]= x3[i-1];
	 
	 rej++;
	 
	 if (rej%nRej == 0) {
	 delta = .9*delta; //VAMOS TER DE ALTERAR ISTO!!!!!!
	 }
	 }
	 
	 else {
	 chi_sq[i] = chi_sq_test;
	 x0[i]= r0l; x1[i]= r1l; x2[i]= r0r; x3[i]= r1r;
	 //rej = 0.; I dont think I can have this here because it is just a propb. acceptance
	 }
	 
	 }
	 
	 outFile << i << " " << chi_sq[i] << " " << x0[i] << " " << x1[i] << " " << x2[i] << " " << x3[i]<< endl;
	 i++;
	 }
	 outFile.close();
	 outFile1.close();
	 
	 double tstart, tstop, ttime;
	 
	 tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
	 
	 tstop = (double)clock()/CLOCKS_PER_SEC;
	 
	 ttime= tstop-tstart; //ttime is how long your code run
	 
	 printf("It took %g time to run the METROPOLIS ALGORITHM for nITER = %d \n", ttime, nIter);
	 
	 */
	
	// MINIMUM CHISQ=ref: (x0,x1,x2,x3) = point_(a,b,c,d)
	
	printf("minimum chiSq  = %g (x0, x1, x2, x3) = (%g, %g, %g,%g)\n", ref, point_a, point_b, point_c, point_d);
	point_a = 0.000103356; point_b= 0.809456; point_c = 916.55; point_d = 701.483;
	r0l= point_a; r1l = point_b;
	r0r = point_c; r1r = point_d;
	bias = funct_elephant_inside_snake(rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
	//writeArrayToFile(rkpc, DoubleTimesArray(100.,bias), name_file+"_bias_final");
	
	mass_true = total_mass_from_bias_elephant_inside_snake(mass_HE, rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
	//writeArrayToFile(rkpc, ArrayDividesDouble(mass_true, oneMsunkg), name_file+"_mass_true_final");
	
	dm_mass = ArrayMinusArray(mass_true, gas_mass);
	//fit to DM density
	fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	//writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_final");
	
	MyFile = name_file + "_chiSQ_r_true_jeans_final";
	outFile.open(MyFile.c_str(),ios::out);
	
	chi_sq_test=0.;
	
	for (int i = 0; i < nBin; i++) {
		double value = pow(mass_true.read(i)-mass_jeans.read(i),2.)/(pow(mass_true.read(i),2.)+ pow(mass_jeans.read(i),2.));
		outFile << rkpc.read(i) << " " << value << endl;
		if (rkpc.read(i) >= 0.001 && rkpc.read(i) <= 10.) {
			chi_sq_test+= value;
		}
	}
	outFile.close();
	printf("minimum chiSq  = %g (x0, x1, x2, x3) = (%g, %g, %g,%g)\n", chi_sq_test, point_a, point_b, point_c, point_d);
	
	
	gas_mass.delete_array(); mass_HE.delete_array(); mass_true.delete_array(); dm_mass.delete_array(); dm_density.delete_array(); mass_jeans.delete_array(); bias.delete_array();
	delete[] fit_parameters;
	//printf("minimum chiSq  = %g (x0, x1, x2, x3) = (%g, %g, %g,%g)\n", ref, point_a, point_b, point_c, point_d);
	
}

double * monte_carlo_loop_hat_min_max_min(string name_file,int nIter, int nRej, double delta, array rkpc, array radius, array gas_density, array gas_mass, array mass_HE, 
										  array disp_velocity, array dm_beta, double * dm_dens_guess, double * x_i, double * min_max_plateau){
	
	ofstream outFile;
	string MyFile;
	int nBin = rkpc.numberBins();
	
	double r0l, r1l, r0r, r1r, min_plateau, max_plateau;
	r0l = x_i[0]; r1l = x_i[1]; r0r = x_i[2]; r1r = x_i[3];
	min_plateau = min_max_plateau[0]; max_plateau = min_max_plateau[1];
	
	array mass_true, dm_mass, dm_density, mass_jeans, bias;
	mass_true.init(nBin); dm_mass.init(nBin); dm_density.init(nBin); mass_jeans.init(nBin); bias.init(nBin);
	
	bias = funct_elephant_inside_snake(rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
	mass_true = total_mass_from_bias_elephant_inside_snake(mass_HE, rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
	dm_mass = ArrayMinusArray(mass_true, gas_mass);
	
	//fit to DM density
	double * fit_parameters; fit_parameters = new double [8];
	fit_parameters = fit_mass_from_density_PL_return_parameters_set_alpha_beta(radius, dm_mass, dm_dens_guess);
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	//	writeArrayToFile(rkpc, ArrayDividesDouble(mass_HE, oneMsunkg), name_file+"_mass_Jeans_w_ansatz");
	
	double chi_sq_test=0.;
	
	for (int i = 0; i < nBin; i++) {
		double value = pow(mass_true.read(i)-mass_jeans.read(i),2.)/(pow(0.1*mass_true.read(i),2.)+ pow(0.1*mass_jeans.read(i),2.));
		//		outFile << rkpc.read(i) << " " << value << endl;
		//if (rkpc.read(i) >= 0.001 && rkpc.read(i) <= 10.) {
		chi_sq_test+= value;
		//}
	}
	
	// INIT METROPOLIS
	
	double chi_sq[nIter], x0[nIter], x1[nIter], x2[nIter], x3[nIter];
	double point_a,point_b, point_c, point_d, ref;
	ref = chi_sq_test;
	chi_sq[0] = chi_sq_test;
	x0[0]= r0l; x1[0]= r1l; x2[0]= r0r; x3[0]= r1r;
	
	MyFile = name_file + "_chiSQ_vs_X0123.dat";
	outFile.open(MyFile.c_str(),ios::out);
	outFile << 0 << " " << chi_sq[0] << " " << x0[0] << " " << x1[0] << " " << x2[0] << " " << x3[0]<< endl;
	
	int i =1; int rej = 0.;
	while (i<nIter){
		int j = 0;
		printf("METROPOLIS: New configuration %d previous chiSQ = %g\n", i, chi_sq[i-1] );
		bool make_new = true;
		while (make_new == true) {
			
			r0l = x0[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			r1l = x1[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			r0r = x2[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			r1r = x3[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			
			// DO we need to set a constraint on r1l and r1r???
			
			mass_true = total_mass_from_bias_elephant_inside_snake(mass_HE, rkpc, r0l, r1l, r0r, r1r, min_plateau, max_plateau);
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (mass_true.read(i) < gas_density.read(i)) make_new = true;
			}
			printf("Rdm number \n");
			j++;
			if (j>30) break;
		}
		if (j>30){
			break;
		}
		dm_mass = ArrayMinusArray(mass_true, gas_density);
		
		fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
		
		dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1],
									 fit_parameters[2],fit_parameters[3]);
		
		mass_jeans= MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
		
		chi_sq_test=0.;
		for (int j = 0; j < nBin; j++) {
			//	if (rkpc.read(j) >= 0.001 && rkpc.read(j) <= 10.) {
			chi_sq_test+= pow(mass_true.read(j)-mass_jeans.read(j),2.)/(pow(mass_true.read(j),2.)+ pow(mass_jeans.read(j),2.));
			//}
		}
		
		if(chi_sq_test < chi_sq[i-1]){ // accepts
			chi_sq[i] = chi_sq_test;
			x0[i]= r0l; x1[i]= r1l; 
			x2[i]= r0r; x3[i]= r1r;
			if (chi_sq_test < ref) {
				point_a = r0l; point_b = r1l;
				point_c = r0r; point_d = r1r;
				ref = chi_sq_test;
			}
			
			rej = 0.;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[i-1]));
			if (r > w) { // i do not accept the new configuration
				chi_sq[i] = chi_sq[i-1];
				x0[i]= x0[i-1]; x1[i]= x1[i-1]; x2[i]= x2[i-1]; x3[i]= x3[i-1];
				
				rej++;
				
				if (rej%nRej == 0) {
					delta = (1.-(1/chi_sq_test))*delta; //VAMOS TER DE ALTERAR ISTO!!!!!!
				}
			}
			
			else {
				chi_sq[i] = chi_sq_test;
				x0[i]= r0l; x1[i]= r1l; x2[i]= r0r; x3[i]= r1r;
				//rej = 0.; I dont think I can have this here because it is just a propb. acceptance
			}
			
		}
		
		outFile << i << " " << chi_sq[i] << " " << x0[i] << " " << x1[i] << " " << x2[i] << " " << x3[i]<< endl;
		i++;
	}
	outFile.close();
	
	// MINIMUM CHISQ=ref: (x0,x1,x2,x3) = point_(a,b,c,d)
	printf("minimum chiSq  = %g (x0, x1, x2, x3) = (%g, %g, %g,%g)\n", ref, point_a, point_b, point_c, point_d);
	
	double * value_param; value_param = new double [5];
	value_param[0] = ref;
	value_param[1] = point_a; value_param[2] = point_b; value_param[3] = point_c; value_param[4] = point_d;
	
	mass_true.delete_array(); dm_mass.delete_array(); dm_density.delete_array(); mass_jeans.delete_array(); bias.delete_array();
	delete[] fit_parameters;
	return value_param;
}

void  genetic_algorithm_same_delta_hat_min_max(string name_file, int nIter_gAl,int nIter, int nRej, double delta, double * x_i, double * min_max_plateau,
											   array rkpc, array radius, array gas_density, array gas_temperature, 
											   array disp_velocity, array dm_beta, double * dm_dens_guess) {
	
	ofstream outFile;
	string MteCarlo_name_file;
	
	// SET INIT - IMAGINE M_HE = M_JEANS. TRUE?
	int nBin = rkpc.numberBins();
	
	array gas_mass, mass_HE; gas_mass.init(nBin); mass_HE.init(nBin);
	gas_mass = massFromDensityAssumeSpherical(radius, gas_density);
	mass_HE = MassFromHE(radius,gas_density,gas_temperature);
	
	double chi_sq_test;
	
	double * values; values = new double [5];
	double * father, * mother;
	father = new double [5]; mother = new double [5];
	
	double * x_i_family, * min_max_family;
	x_i_family = new double [4]; min_max_family = new double [2];
	
	for (int i = 0; i < 4; i++) x_i_family[i] = x_i[i];
	MteCarlo_name_file = name_file + "_F1";
	father = monte_carlo_loop_hat_min_max_min(MteCarlo_name_file,nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
											  disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_plateau);
	
	mother[0] = 180.;
	for (int i = 1; i < 5; i++) mother[i] = x_i[i];
	
	int iter = 0;
	while (iter < nIter_gAl) {
		char s[0];
		sprintf(s, "%d", iter );
		double * son1, * son2, * son3, * son4;
		son1 = new double [5];  son2 = new double [5]; son3 = new double [5]; son4 = new double [5];
		
		for (int i = 0; i < 2; i++) x_i_family[i] = father[i+1];
		for (int i = 2; i < 4; i++) x_i_family[i] = mother[i+1];
		min_max_family[0]  = min_max_plateau[0]; min_max_family[1] = min_max_plateau[1];
		MteCarlo_name_file = name_file + "_S1_" + s;
		son1 = monte_carlo_loop_hat_min_max_min(MteCarlo_name_file, nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
		printf("Iter = %d, Son1: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son1[0], son1[1], son1[2]);
		
		for (int i = 0; i < 2; i++) x_i_family[i] = mother[i+1];
		for (int i = 2; i < 4; i++) x_i_family[i] = father[i+1];
		min_max_family[0]  = min_max_plateau[0]; min_max_family[1] = min_max_plateau[1];
		MteCarlo_name_file = name_file + "_S2_" + s;
		son2 = monte_carlo_loop_hat_min_max_min(MteCarlo_name_file, nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
		printf("Iter = %d, Son2: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son2[0], son2[1], son2[2]);
		
		for (int i = 0; i < 2; i++) x_i_family[i] = mother[i+2];
		for (int i = 2; i < 4; i++) x_i_family[i] = father[i+1];
		min_max_family[0]  = 0.5 * min_max_plateau[0]; min_max_family[1] = 0.5 * min_max_plateau[1];
		MteCarlo_name_file = name_file + "_S3_" + s;
		
		son3 = monte_carlo_loop_hat_min_max_min(MteCarlo_name_file,nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
		printf("Iter = %d, Son3: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son3[0], son3[1], son3[2]);
		
		for (int i = 0; i < 2; i++) x_i_family[i] = mother[i+1];
		for (int i = 2; i < 4; i++) x_i_family[i] = father[i-1];
		min_max_family[0]  = 0.5 * min_max_plateau[0]; min_max_family[1] = 0.5 * min_max_plateau[1];
		MteCarlo_name_file = name_file + "_S4_" + s;
		
		son4 = monte_carlo_loop_hat_min_max_min(MteCarlo_name_file,nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
		printf("Iter = %d, Son4: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son4[0], son4[1], son4[2]);
		
		if (son1[0] < father[0]) {
			for (int i = 0; i < 5; i++) father[i] = son1[i];
		}
		if (son1[0] > father[0] && son1[0] < mother[0]) {
			for (int i = 0; i < 5; i++) mother[i] = son1[i];
		}
		if (son2[0] < father[0]) {
			for (int i = 0; i < 5; i++) father[i] = son2[i];
		}
		if (son2[0] > father[0] && son2[0] < mother[0]) {
			for (int i = 0; i < 5; i++) mother[i] = son2[i];
		}
		if (son3[0] < father[0]) {
			for (int i = 0; i < 5; i++) father[i] = son3[i];
		}
		if (son3[0] > father[0] && son3[0] < mother[0]) {
			for (int i = 0; i < 5; i++) mother[i] = son3[i];
		}
		if (son4[0] < father[0]) {
			for (int i = 0; i < 5; i++) father[i] = son4[i];
		}
		if (son4[0] > father[0] && son4[0] < mother[0]) {
			for (int i = 0; i < 5; i++) mother[i] = son4[i];
		}
		
		printf("END Father: Best ChiSQ = %g for (x0,x1)=(%g,%g) (x2,x3)=(%g,%g)\n", father[0], father[1], father[2], father[3], father[4]);
		printf("END Mother: Best ChiSQ = %g for (x0,x1)=(%g,%g) (x2,x3)=(%g,%g)\n", mother[0], mother[1], mother[2],  mother[3], mother[4]);
		
		iter++;
		delete[] son1; delete[] son2; delete[] son3; delete[] son4;
	}
	
	printf("END Father: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", father[0], father[1], father[2]);
	printf("END Mother: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", mother[0], mother[1], mother[2]);
	
	gas_mass.delete_array(); mass_HE.delete_array();
	delete[] values;
	delete[] x_i_family; delete[] min_max_family;
	
	
}

double * monte_carlo_loop_step_min_max(string name_file,int nIter, int nRej, double delta, array rkpc, array radius, array gas_density, array gas_mass, array mass_HE, 
									   array disp_velocity, array dm_beta, double * dm_dens_guess, double * x_i, double * min_max_plateau){
	
	ofstream outFile, outFile1;
	string MyFile;
	int nBin = rkpc.numberBins();
	
	double r0l, r1l, r0r, r1r, min_plateau, max_plateau;
	r0l = x_i[0]; r1l = x_i[1];
	min_plateau = min_max_plateau[0]; max_plateau = min_max_plateau[1];
	
	array mass_true, dm_mass, dm_density, mass_jeans, bias;
	mass_true.init(nBin); dm_mass.init(nBin); dm_density.init(nBin); mass_jeans.init(nBin); bias.init(nBin);
	
	dm_mass = ArrayMinusArray(mass_HE, gas_mass);
	
	//fit to DM density
	double * fit_parameters; fit_parameters = new double [8];
	//fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	fit_parameters = fit_mass_from_density_PL_return_parameters_set_alpha_beta(radius, dm_mass, dm_dens_guess);
	
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	writeArrayToFile(radius, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_init");
	
	
	//bias = funct_min_to_max(rkpc, r0l, r1l, min_plateau, max_plateau);
	//mass_true = total_mass_from_bias(mass_HE, bias);
	//dm_mass = ArrayMinusArray(mass_true, gas_mass);
	
	//fit to DM density
	//double * fit_parameters; fit_parameters = new double [8];
	//fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	//fit_parameters = fit_mass_from_density_PL_return_parameters_set_alpha_beta(radius, dm_mass, dm_dens_guess);
	
	//dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	//mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	//writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_ans");
	double delta_ref = delta;
	
	double chi_sq_test=0.;
	int k =0;
	for (int i = 0; i < nBin; i++) {
		//double value = pow(mass_true.read(i)-mass_jeans.read(i),2.)/(pow(mass_true.read(i),2.)+ pow(mass_jeans.read(i),2.));
		double value = pow(mass_true.read(i)-mass_jeans.read(i),2.)/(pow(0.1*mass_true.read(i),2.)+ pow(0.1*mass_jeans.read(i),2.));
		//		outFile << rkpc.read(i) << " " << value << endl;
		if (rkpc.read(i) >= 0.001 && rkpc.read(i) <= 10.) {
			chi_sq_test+= value;
			k++;
		}
	}
	chi_sq_test = chi_sq_test/k;
	printf("chiSQ init = %g \n", chi_sq_test);
	// INIT METROPOLIS
	
	double chi_sq[nIter], x0[nIter], x1[nIter], x2[nIter], x3[nIter];
	double point_a, point_b, ref;
	ref = chi_sq_test;
	chi_sq[0] = chi_sq_test;
	x0[0]= r0l; x1[0]= r1l;
	
	MyFile = name_file + "_best_chiSQ_vs_X0.dat";
	outFile1.open(MyFile.c_str(),ios::out);
	
	double chi_rej = chi_sq[0];
	double tstart, tstop, ttime;
	
	int i =1; int rej = 0.; int accept = 0;
	char num_string[nIter];
	while (i<nIter){
		sprintf(num_string, "%d", i);
		MyFile = name_file + "_"+num_string;
		
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		
		int j = 0;
		printf("METROPOLIS: New configuration %d previous chiSQ = %g when (x0,x1) = (%g,%g) delta = %g\n", i, chi_sq[i-1], x0[i-1], x1[i-1], delta);
		bool make_new = true;
		while (make_new == true) {
			
			r0l = x0[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			r1l = x1[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			
			printf("r0l = %g and r1l = %g \n", r0l, r1l);
			bias = funct_min_to_max(rkpc, r0l, r1l, min_plateau, max_plateau);
			mass_true = total_mass_from_bias(mass_HE, bias);
			
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (mass_true.read(i) < gas_density.read(i)) make_new = true;
			}
			printf("Rdm number \n");
			j++;
			if (j>30) break;
		}
		if (j>30){
			break;
		}
		
		dm_mass = ArrayMinusArray(mass_true, gas_density);
		
		fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
		
		dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1],
									 fit_parameters[2],fit_parameters[3]);
		mass_jeans= MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
		
		chi_sq_test=0.;
		k=0;
		for (int j = 0; j < nBin; j++) {
			if (rkpc.read(j) >= 0.001 && rkpc.read(j) <= 10.) {
				chi_sq_test+= pow(mass_true.read(j)-mass_jeans.read(j),2.)/(pow(0.1*mass_true.read(j),2.)+ pow(0.1*mass_jeans.read(j),2.));
				//printf("calculate chiSQ\n");
				k++;
			}
		}
		chi_sq_test = chi_sq_test/k;
		
		if(chi_sq_test < chi_sq[i-1]){ // accepts
			chi_sq[i] = chi_sq_test;
			x0[i]= r0l;
			x1[i]= r1l;
			//writeArrayToFile(radius, ArrayDividesDouble(mass_true,oneMsunkg), MyFile + "_true");
			
			//writeArrayToFile(radius, ArrayDividesDouble(mass_jeans,oneMsunkg), MyFile+ "_jeans");
			printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[i-1]);
			if (chi_sq_test < ref) {
				point_a = r0l;
				point_b = r1l;
				ref = chi_sq_test;
				outFile1 << i << " " << ref << " " << point_a << " " << point_b << endl;
			}
			//rej = 0.;
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[i-1]));
			printf("r = %g   w = %g \n", r, w);
			if (r > w) { // i do not accept the new configuration
				chi_sq[i] = chi_sq[i-1];
				x0[i]= x0[i-1]; 
				x1[i]= x1[i-1]; 
				
				rej++;
				
				if (rej%nRej == 0) {
					
					delta = 0.85*delta; //VAMOS TER DE ALTERAR ISTO!!!!!!
					chi_rej = chi_sq_test;
				}
			}
			
			if (r<w) {
				chi_sq[i] = chi_sq_test;
				x0[i]= r0l; 
				x1[i]= r1l; 
				accept++;
			}
			
		}
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		printf("It took [%g] time to run one LOOP. nIter = %d  nAccept = %d\n", ttime, nIter, accept);
		
		i++;
	}
	
	MyFile = name_file + "_chiSQ_vs_X01.dat";
	outFile.open(MyFile.c_str(),ios::out);
	
	for (int i = 0; i < nIter; i++) {
		outFile << i << " " << chi_sq[i] << " " << x0[i] << " " << x1[i] << endl;
		
	}
	outFile.close();
	
	
	mass_true.delete_array(); dm_mass.delete_array(); dm_density.delete_array(); mass_jeans.delete_array();bias.delete_array();
	
	bias = funct_min_to_max(rkpc, point_a, point_b, min_plateau, max_plateau);
	mass_true = total_mass_from_bias(mass_HE, bias);
	dm_mass = ArrayMinusArray(mass_true, gas_mass);
	writeArrayToFile(rkpc, ArrayDividesDouble(mass_true, oneMsunkg), name_file+"_mass_true_final");
	//fit to DM density
	//fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_final_free");
	
	fit_parameters = fit_mass_from_density_PL_return_parameters_set_alpha_beta(radius, dm_mass, dm_dens_guess);
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	writeArrayToFile(rkpc, Gamma(radius, dm_density), name_file+"_gamma_dm_final");
	writeArrayToFile(rkpc, ArrayDividesDouble(dm_density,atom), name_file+"_dm_dens_final");
	
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_final_nfw");
	writeArrayToFile(rkpc, Gamma(radius, dm_density), name_file+"_gamma_dm_nfw_final");
	writeArrayToFile(rkpc, ArrayDividesDouble(dm_density,atom), name_file+"_dm_dens_nfw_final");
	
	delete [] fit_parameters;
	
	// MINIMUM CHISQ=ref: (x0,x1,x2,x3) = point_(a,b,c,d)
	printf("minimum chiSq  = %g (x0, x1) = (%g, %g)\n", ref, point_a, point_b);
	
	double * value_param; value_param = new double [3];
	value_param[0] = ref;
	value_param[1] = point_a; value_param[2] = point_b;
	
	return value_param;
}


double * monte_carlo_loop_step_max_min(string name_file,int nIter, int nRej, double delta, array rkpc, array radius, array gas_density, array gas_mass, array mass_HE, 
									   array disp_velocity, array dm_beta, double * dm_dens_guess, double * x_i, double * min_max_plateau){
	
	ofstream outFile, outFile1;
	string MyFile;
	int nBin = rkpc.numberBins();
	
	double r0l, r1l, r0r, r1r, min_plateau, max_plateau;
	r0l = x_i[0]; r1l = x_i[1];
	min_plateau = min_max_plateau[0]; max_plateau = min_max_plateau[1];
	
	array mass_true, dm_mass, dm_density, mass_jeans, bias;
	mass_true.init(nBin); dm_mass.init(nBin); dm_density.init(nBin); mass_jeans.init(nBin); bias.init(nBin);
	
	dm_mass = ArrayMinusArray(mass_HE, gas_mass);
	
	//fit to DM density
	double * fit_parameters; fit_parameters = new double [8];
	//fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	fit_parameters = fit_mass_from_density_PL_return_parameters_set_alpha_beta(radius, dm_mass, dm_dens_guess);
	
	//dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	//mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	//writeArrayToFile(radius, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_init");
	
	
	//bias = funct_min_to_max(rkpc, r0l, r1l, min_plateau, max_plateau);
	//mass_true = total_mass_from_bias(mass_HE, bias);
	//dm_mass = ArrayMinusArray(mass_true, gas_mass);
	
	//fit to DM density
	//double * fit_parameters; fit_parameters = new double [8];
	//fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	//fit_parameters = fit_mass_from_density_PL_return_parameters_set_alpha_beta(radius, dm_mass, dm_dens_guess);
	
	//dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	//mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	//writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_ans");
	double delta_ref = delta;
	
	double chi_sq_test=0.;
	int k =0;
	//for (int i = 0; i < nBin; i++) {
	//double value = pow(mass_true.read(i)-mass_jeans.read(i),2.)/(pow(mass_true.read(i),2.)+ pow(mass_jeans.read(i),2.));
	//	double value = pow(mass_true.read(i)-mass_jeans.read(i),2.)/(pow(0.1*mass_true.read(i),2.)+ pow(0.1*mass_jeans.read(i),2.));
	//		outFile << rkpc.read(i) << " " << value << endl;
	//	if (rkpc.read(i) >= 0.001 && rkpc.read(i) <= 10.) {
	//		chi_sq_test+= value;
	//		k++;
	//	}
	//} 
	k=nBin;
	chi_sq_test=nBin;
	chi_sq_test = chi_sq_test/k;
	printf("chiSQ init = %g \n", chi_sq_test);
	// INIT METROPOLIS
	
	double chi_sq[nIter], x0[nIter], x1[nIter], x2[nIter], x3[nIter];
	double point_a, point_b, ref;
	ref = chi_sq_test;
	chi_sq[0] = chi_sq_test;
	x0[0]= r0l; x1[0]= r1l;
	
	MyFile = name_file + "_best_chiSQ_vs_X0.dat";
	outFile1.open(MyFile.c_str(),ios::out);
	
	double chi_rej = chi_sq[0];
	double tstart, tstop, ttime;
	
	int i =1; int rej = 0.; int accept = 0;
	char num_string[nIter];
	while (i<nIter){
		sprintf(num_string, "%d", i);
		MyFile = name_file + "_"+num_string;
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		
		int j = 0;
		printf("METROPOLIS: New configuration %d previous chiSQ = %g when (x0,x1) = (%g,%g) delta = %g\n", i, chi_sq[i-1], x0[i-1], x1[i-1], delta);
		bool make_new = true;
		while (make_new == true) {
			
			r0l = x0[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			r1l = x1[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			
			printf("r0l = %g and r1l = %g \n", r0l, r1l);
			bias = funct_max_to_min(rkpc, r0l, r1l, min_plateau, max_plateau);
			mass_true = total_mass_from_bias(mass_HE, bias);
			
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (mass_true.read(i) < gas_density.read(i)) make_new = true;
			}
			printf("Rdm number \n");
			j++;
			if (j>30) break;
		}
		if (j>30){
			break;
		}
		
		dm_mass = ArrayMinusArray(mass_true, gas_density);
		
		//fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
		
		//dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1],
		//							 fit_parameters[2],fit_parameters[3]);
		dm_density = DensityFromMassSpherical(radius, dm_mass);
		mass_jeans= MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
		
		chi_sq_test=0.;
		k=0;
		for (int j = 0; j < nBin; j++) {
			if (rkpc.read(j) >= 0.001 && rkpc.read(j) <= 10.) {
				chi_sq_test+= pow(mass_true.read(j)-mass_jeans.read(j),2.)/(pow(0.1*mass_true.read(j),2.)+ pow(0.1*mass_jeans.read(j),2.));
				//printf("calculate chiSQ\n");
				k++;
			}
		}
		chi_sq_test = chi_sq_test/k;
		
		if(chi_sq_test < chi_sq[i-1]){ // accepts
			chi_sq[i] = chi_sq_test;
			x0[i]= r0l;
			x1[i]= r1l;
			//writeArrayToFile(radius, ArrayDividesDouble(mass_true,oneMsunkg), MyFile + "_true");
			
			//writeArrayToFile(radius, ArrayDividesDouble(mass_jeans,oneMsunkg), MyFile+ "_jeans");
			printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[i-1]);
			if (chi_sq_test < ref) {
				point_a = r0l;
				point_b = r1l;
				ref = chi_sq_test;
				outFile1 << i << " " << ref << " " << point_a << " " << point_b << endl;
			}
			//rej = 0.;
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[i-1]));
			printf("r = %g   w = %g \n", r, w);
			if (r > w) { // i do not accept the new configuration
				chi_sq[i] = chi_sq[i-1];
				x0[i]= x0[i-1]; 
				x1[i]= x1[i-1]; 
				
				rej++;
				
				if (rej%nRej == 0) {
					
					delta = 0.85*delta; //VAMOS TER DE ALTERAR ISTO!!!!!!
					chi_rej = chi_sq_test;
				}
			}
			
			if (r<w) {
				chi_sq[i] = chi_sq_test;
				x0[i]= r0l; 
				x1[i]= r1l; 
				accept++;
			}
			
		}
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		printf("It took [%g] time to run one LOOP. nIter = %d  nAccept = %d\n", ttime, nIter, accept);
		
		i++;
	}
	
	MyFile = name_file + "_chiSQ_vs_X01.dat";
	outFile.open(MyFile.c_str(),ios::out);
	
	for (int i = 0; i < nIter; i++) {
		outFile << i << " " << chi_sq[i] << " " << x0[i] << " " << x1[i] << endl;
		
	}
	outFile.close();
	
	
	mass_true.delete_array(); dm_mass.delete_array(); dm_density.delete_array(); mass_jeans.delete_array();bias.delete_array();
	
	bias = funct_min_to_max(rkpc, point_a, point_b, min_plateau, max_plateau);
	mass_true = total_mass_from_bias(mass_HE, bias);
	dm_mass = ArrayMinusArray(mass_true, gas_mass);
	writeArrayToFile(rkpc, ArrayDividesDouble(mass_true, oneMsunkg), name_file+"_mass_true_final");
	//fit to DM density
	//fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	fit_parameters = fit_mass_from_density_PL_return_parameters(radius, dm_mass, dm_dens_guess);
	
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_final_free");
	
	fit_parameters = fit_mass_from_density_PL_return_parameters_set_alpha_beta(radius, dm_mass, dm_dens_guess);
	dm_density = DensityPowerLaw(radius, fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	writeArrayToFile(rkpc, Gamma(radius, dm_density), name_file+"_gamma_dm_final");
	writeArrayToFile(rkpc, ArrayDividesDouble(dm_density,atom), name_file+"_dm_dens_final");
	
	mass_jeans = MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	writeArrayToFile(rkpc, ArrayDividesDouble(mass_jeans, oneMsunkg), name_file+"_mass_Jeans_final_nfw");
	writeArrayToFile(rkpc, Gamma(radius, dm_density), name_file+"_gamma_dm_nfw_final");
	writeArrayToFile(rkpc, ArrayDividesDouble(dm_density,atom), name_file+"_dm_dens_nfw_final");
	
	delete [] fit_parameters;
	
	// MINIMUM CHISQ=ref: (x0,x1,x2,x3) = point_(a,b,c,d)
	printf("minimum chiSq  = %g (x0, x1) = (%g, %g)\n", ref, point_a, point_b);
	
	double * value_param; value_param = new double [3];
	value_param[0] = ref;
	value_param[1] = point_a; value_param[2] = point_b;
	
	return value_param;
}



void  genetic_algorithm_same_delta_step_min_max(string name_file, int nIter_gAl,int nIter, int nRej, double delta, double * x_i, double * min_max_plateau,
												array rkpc, array radius, array gas_density, array gas_temperature, 
												array disp_velocity, array dm_beta, double * dm_dens_guess) {
	
	ofstream outFile;
	string MteCarlo_name_file;
	
	// SET INIT - IMAGINE M_HE = M_JEANS. TRUE?
	int nBin = rkpc.numberBins();
	
	array gas_mass, mass_HE; gas_mass.init(nBin); mass_HE.init(nBin);
	gas_mass = massFromDensityAssumeSpherical(radius, gas_density);
	mass_HE = MassFromHE(radius,gas_density,gas_temperature);
	
	double chi_sq_test;
	
	double * values; values = new double [3];
	double * father, * mother;
	
	father = new double [4]; mother = new double [4];
	
	double * x_i_family, * min_max_family;
	x_i_family = new double [2]; min_max_family = new double [2];
	
	for (int i = 0; i < 2; i++) x_i_family[i] = x_i[i];
	MteCarlo_name_file = name_file + "_F1";
	father = monte_carlo_loop_step_min_max(MteCarlo_name_file, nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
										   disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_plateau);
	
	mother[0] = 180.;
	for (int i = 1; i < 3; i++) mother[i] = x_i[i];
	mother [3] = delta;
	
	
	
	int iter = 0;
	while (iter < nIter_gAl) {
		double * son1, * son2, * son3, * son4;
		son1 = new double [4];  son2 = new double [4]; son3 = new double [4]; son4 = new double [4];
		
		bool S1, S2, S3, S4;
		S1=S2=S3=S4=false;
		char s[0];
		sprintf(s, "%d", iter );
		if (father[1] < mother[2]) {
			x_i_family[0] = father[1];
			x_i_family[1] = mother[2];
			min_max_family[0]  = min_max_plateau[0]; min_max_family[1] = min_max_plateau[1];
			
			MteCarlo_name_file = name_file + "_S1_" + s;
			son1 = monte_carlo_loop_step_min_max(MteCarlo_name_file, nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												 disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
			printf("Iter = %d, Son1: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son1[0], son1[1], son1[2]);
			S1 = true;
		}
		if (mother[1] < father[2]) {
			x_i_family[0] = mother[1];
			x_i_family[1] = father[2];
			min_max_family[0]  = min_max_plateau[0]; min_max_family[1] = min_max_plateau[1];
			
			MteCarlo_name_file = name_file + "_S2_" + s;
			son2 = monte_carlo_loop_step_min_max(MteCarlo_name_file, nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												 disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
			printf("Iter = %d, Son2: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son2[0], son2[1], son2[2]);
			
			S2 = true;
		}
		
		if (mother[2] < father[1]) {
			// if mother(2) > father(1) this will be consider a mutation it is not even a function
			x_i_family[0] = mother[2];
			x_i_family[1] = father[1];
			min_max_family[0]  = 0.5 * min_max_plateau[0]; min_max_family[1] = 0.5 * min_max_plateau[1];
			
			MteCarlo_name_file = name_file + "_S3_" + s;
			son3 = monte_carlo_loop_step_min_max(MteCarlo_name_file, nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												 disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
			printf("Iter = %d, Son3: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son3[0], son3[1], son3[2]);
			
			S3 = true;
		}
		
		if (father [2] < mother[1]) {
			
			MteCarlo_name_file = name_file + "_S4_" + s;
			x_i_family[0] = father[2];
			x_i_family[1] = mother[1];
			min_max_family[0]  = 0.5 * min_max_plateau[0]; min_max_family[1] = 0.5 * min_max_plateau[1];
			son4 = monte_carlo_loop_step_min_max(MteCarlo_name_file, nIter, nRej, delta, rkpc, radius, gas_density, gas_temperature, mass_HE,
												 disp_velocity, dm_beta, dm_dens_guess, x_i_family, min_max_family);
			
			printf("Iter = %d, Son3: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", iter, son4[0], son4[1], son4[2]);
			S4 = true;
		}
		
		if (son1[0] < father[0] && S1 == true) {
			for (int i = 0; i < 4; i++) father[i] = son1[i];
		}
		if (son1[0] > father[0] && son1[0] < mother[0]  && S1 == true) {
			for (int i = 0; i < 3; i++) mother[i] = son1[i]; //delta for the mother will always start at the previous, 
			//we are only doing this for the father to see if we get faster to the minima
		}
		if (son2[0] < father[0]  && S2 == true) {
			for (int i = 0; i < 4; i++) father[i] = son2[i];
		}
		if (son2[0] > father[0] && son2[0] < mother[0]  && S2 == true) {
			for (int i = 0; i < 3; i++) mother[i] = son2[i];
		}
		if (son3[0] < father[0]  && S3 == true) {
			for (int i = 0; i < 4; i++) father[i] = son3[i];
		}
		if (son3[0] > father[0] && son3[0] < mother[0]  && S3 == true) {
			for (int i = 0; i < 3; i++) mother[i] = son3[i];
		}
		if (son4[0] < father[0]  && S4 == true) {
			for (int i = 0; i < 4; i++) father[i] = son4[i];
		}
		if (son4[0] > father[0] && son4[0] < mother[0]  && S4 == true) {
			for (int i = 0; i < 3; i++) mother[i] = son4[i];
		}
		delete []son1;
		delete []son2;
		delete []son3;
		delete []son4;
		
		
		iter++;
	}
	
	printf("END Father: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", father[0], father[1], father[2]);
	printf("END Mother: Best ChiSQ = %g for (x0,x1)=(%g,%g)\n", mother[0], mother[1], mother[2]);
	
}

array * Metropolis_get_jeans_and_delta_m(string name_file, int u, int nIter, int nRej, double delta, array radius, array gas_dens, array gas_mass, array he_mass,
										 array disp_velocity, array dm_beta, 
										 double * dm_dens_guess_fit, double * x){
	
	printf("Starting Metropolis\n");
	int nBin = radius.numberBins();
	
	array * jeans_delta_m; jeans_delta_m = new array [2];
	jeans_delta_m[0].init(nBin); jeans_delta_m[1].init(nBin);
	
	double mass_true_testing[nBin];
	array test_true_mass; test_true_mass.init(nBin);
	
	double chi_sq[nIter]; chi_sq[0] = 100.; int accept = 0; int rej = 0;
	
	array dm_mass, dm_density, mass_jeans, bias;
	dm_mass.init(nBin); dm_density.init(nBin); mass_jeans.init(nBin); bias.init(nBin);
	
	array ref_true_mass; ref_true_mass.init(nBin);
	
	char s[1];
	sprintf(s, "%d", u );
	
	ofstream outFile, outFile1;
	string MyFile;
	MyFile = name_file+s+"_normal_chiSQ_vs_points.dat";
	outFile.open(MyFile.c_str(),ios::out);
	MyFile =  name_file+s+"_best_chiSQ_vs_points.dat";
	outFile1.open(MyFile.c_str(),ios::out);
	
	double tstart, tstop, ttime;
	
	//double * fit_parameters; fit_parameters = new double [8];
	
	double  x0[nIter], x1[nIter], x2[nIter], x3[nIter];
	double point_a, point_b, point_c, point_d, ref;
	ref = chi_sq[0];
	x0[0]= x[0]; x1[0]= x[1];x2[0]= x[2]; x3[0]= x[3];
	
	double xx0, xx1, xx2, xx3;
	
	double min_plateau = -.9;
	double max_plateau = 0.;
	
	int i = 1; int k;
	while (i< nIter) {
		printf("Starting iteration N = %d\n", i);
		tstart = (double)clock()/CLOCKS_PER_SEC; 
		bool make_new = false;
		k=0;
		while (make_new == false) {
			k++;
			if (k>30) break;
			for (int j = 0; j < nBin; j++) {
				xx0 =  x0[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
				xx1 =  x1[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
				xx2 =  x2[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
				xx3 =  x3[i-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
				//test_true_mass.set(j, true_mass.read(j) *  pow(10.,cat_rng_uniform_double_min_max(-delta,delta)));
				
				printf("(x0, x1, x2, x3) = (%g, %g, %g, %g) \n", xx0, xx1, xx2,xx3);
				bias = funct_elephant_inside_snake(radius, xx0, xx1, xx2, xx3, min_plateau, max_plateau);
				test_true_mass = total_mass_from_bias(he_mass, bias);
			}
			make_new = true;
			for (int j = 0; j < nBin-1; j++) {
				if (test_true_mass.read(j+1)<=test_true_mass.read(j)) {
					make_new = false;
					printf("a mass nao cresce i = %d\n", j);
				}
			}
			for (int j =0 ; j < nBin; j++) {
				if (test_true_mass.read(j)<=gas_mass.read(j)) {
					make_new = false;
					printf("a massa total e menos que a mass do gas i = %d\n",j);
				}
			}
			dm_mass = ArrayMinusArray(test_true_mass, gas_mass);
			dm_density = DensityFromMassSpherical(radius, dm_mass);
			
			for (int j =0 ; j < nBin-1; j++) {
				if (dm_density.read(j)<=dm_density.read(j+1)) {
					make_new = false;
					printf("a dens dm esta a subir i = %d\n",j);
				}
			}
			
			
		}
		if (k > 30) break;
		printf("Has found a total mass to test\n");
		mass_jeans= MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
		
		/*		writeArrayToFile(radius, test_true_mass, "_test_true_mass_Sat");
		 writeArrayToFile(radius, he_mass, "_he_mass_Sat");
		 writeArrayToFile(radius, gas_mass, "_gas_mass_Sat");
		 writeArrayToFile(radius, dm_mass, "_dm_mass_Sat");		
		 writeArrayToFile(radius, dm_density, "_dm_dens_Sat");
		 writeArrayToFile(radius, mass_jeans, "_jeans_mass_Sat");
		 writeArrayToFile(radius, ArraySQ(disp_velocity), "_disp_vel_Sat");
		 //*/		
		//cin >> make_new;
		
		
		double chi_sq_test=0.;
		for (int j = 0; j < nBin; j++) {
			chi_sq_test+= pow(test_true_mass.read(j)-mass_jeans.read(j),2.)/(pow(0.1*test_true_mass.read(j),2.)+ pow(0.1*mass_jeans.read(j),2.));
		}
		
		chi_sq_test = chi_sq_test/nBin;
		printf("chi_sq_test = %g\n", chi_sq_test);
		
		if(chi_sq_test < chi_sq[i-1]){ // accepts
			chi_sq[i] = chi_sq_test;
			x0[i]= xx0;
			x1[i]= xx1;
			x2[i]= xx2;
			x3[i]= xx3;
			
			//writeArrayToFile(radius, ArrayDividesDouble(mass_true,oneMsunkg), MyFile + "_true");
			
			//writeArrayToFile(radius, ArrayDividesDouble(mass_jeans,oneMsunkg), MyFile+ "_jeans");
			printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[i-1]);
			if (chi_sq_test < ref) {
				point_a = xx0;
				point_b = xx1;
				point_c = xx2;
				point_d = xx3;
				
				ref = chi_sq_test;
				outFile1 << i << " " << ref << " " << point_a << " " << point_b << " " << point_c<<" " << point_d<< endl;
			}
			//rej = 0.;
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[i-1]));
			printf("r = %g   w = %g \n", r, w);
			if (r > w) { // i do not accept the new configuration
				chi_sq[i] = chi_sq[i-1];
				x0[i]= x0[i-1]; 
				x1[i]= x1[i-1]; 
				x2[i]= x2[i-1]; 
				x3[i]= x3[i-1]; 
				
				rej++;
				
				if (rej%nRej == 0) {
					
					delta = 0.85*delta; //VAMOS TER DE ALTERAR ISTO!!!!!!
					//chi_rej = chi_sq_test;
				}
			}
			
			if (r<w) {
				chi_sq[i] = chi_sq_test;
				x0[i]= xx0; 
				x1[i]= xx1; 
				x2[i]= xx2; 
				x3[i]= xx3; 
				accept++;
			}
			
		}
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		printf("It took [%g] time to run one LOOP. nIter = %d  nAccept = %d\n", ttime, nIter, accept);
		
		i++;
	}
	
	outFile.close();
	outFile1.close(); 
	
	printf("chi_sq = %g <<>>>(x0, x1, x2, x3) = (%g, %g, %g, %g) \n", ref, point_a,point_b, point_c, point_d);
	bias = funct_elephant_inside_snake(radius, xx0, xx1, xx2, xx3, min_plateau, max_plateau);
	test_true_mass = total_mass_from_bias(he_mass, bias);
	
	
	
	
	jeans_delta_m[1] = CopyArray(test_true_mass);
	writeArrayToFile(radius, jeans_delta_m[1], "_true_total_mass");
	
	dm_mass = ArrayMinusArray(jeans_delta_m[1], gas_mass);
	dm_density = DensityFromMassSpherical(radius, dm_mass);
	mass_jeans= MassFromJeans(radius, dm_density, ArraySQ(disp_velocity), dm_beta);
	writeArrayToFile(radius, jeans_delta_m[1], "_new_jeans_mass");
	
	writeArrayToFile(radius, ArrayDividesArray(ArrayMinusArray(he_mass,ref_true_mass),ref_true_mass), "_new_jeans_mass");
	
	
	jeans_delta_m[0] = CopyArray(mass_jeans);
	
	
	
	ref_true_mass.delete_array(); mass_jeans.delete_array(); dm_density.delete_array(); dm_mass.delete_array();
	//true_mass.delete_array();
	//mass_true_testing.delete_array();
	
	return jeans_delta_m;
}


double * Metropolis_get_jeans_equal_to_he(array radius, array he_mass, array he_error, array jeans_mass, array jeans_error,
										  array gas_mass, array dm_beta, array dm_disp_vel, double * x_guess, int nIter, int nRej,double delta, string name_file){
	
	int nBin = radius.numberBins();
	
	double * x_i_return; x_i_return = new double [4];
	
	double chi_sq_test_r =0.;
	for (int i = 0; i < nBin; i ++) {
		double value = pow(he_mass.read(i)-jeans_mass.read(i),2.)/(pow(0.1*he_mass.read(i),2)+pow(0.1*jeans_mass.read(i),2.));
		chi_sq_test_r+= value;
	}
	chi_sq_test_r = chi_sq_test_r/nBin;
	printf("INIT >>>>>> take Jeans from HE <<<<<<<< chiSQ = %g\n", chi_sq_test_r);
	int testing;
	//cin >> testing;
	double chi_sq[nIter]; chi_sq[0] = chi_sq_test_r;
	double chi_sq_ref = chi_sq_test_r;
	double point_a, point_b, point_c, point_d;
	
	double x_0[nIter], x_1[nIter], x_2[nIter], x_3[nIter];
	double x0, x1, x2,x3;
	x0 = x_guess[0]; x1 = x_guess[1]; x2 = x_guess[2]; x3=x_guess[3];
	x_0[0] = x0; x_1[0]=x1; x_2[0]=x2; x_3[0]=x3;
	
	double min_plateau = x_guess[4];
	double max_plateau = x_guess[5];
	
	array bias, total_true_mass, dm_density;
	bias.init(nBin); total_true_mass.init(nBin); dm_density.init(nBin); 
	
	double tstart, tstop, ttime;
	
	string MyFile = name_file + "_chi_sq_vs_x01.txt";
	ofstream outFile0; outFile0.open(MyFile.c_str(), ios::out);
	
	MyFile = name_file + "_chi_sq_vs_x23.txt";
	ofstream outFile1; outFile1.open(MyFile.c_str(), ios::out);
	
	
	int iter = 1; int rej = 0.; int accept = 0;
	int sum_rej = 0;
	srand(time(NULL));
	/*while (iter < nIter){
	 printf("New Run\n");
	 printf("New chi_sq = %g, Old chi_sq = %g \n", chi_sq_test, chi_sq[iter-1]);
	 
	 tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
	 
	 bool make_new = true;
	 int sum_rej = 0;
	 while (make_new == true) {
	 x0 = x_0[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 x1 = x_1[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 x2 = x_2[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 x3 = x_3[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
	 
	 bias = funct_elephant_inside_snake(radius, x0, x1, x2, x3, min_plateau, max_plateau);
	 total_true_mass = total_mass_from_bias(he_mass, bias);
	 make_new = false;
	 for (int i = 0; i < nBin; i++) {
	 if (total_true_mass.read(i) < gas_mass.read(i)){
	 make_new = true;
	 //sum_rej++;
	 }
	 }
	 //if(sum_rej>nIter) break;
	 
	 dm_density = DensityFromMassSpherical(radius, ArrayMinusArray (total_true_mass, gas_mass));
	 
	 for (int i = 0; i < nBin-1; i++) {
	 if (dm_density.read(i) < dm_density.read(i+1)){
	 make_new = true;
	 //sum_rej++;
	 }
	 }
	 
	 jeans_mass = MassFromJeans(radius, dm_density, ArraySQ(dm_disp_vel), dm_beta);
	 
	 
	 
	 chi_sq_test = 0.;
	 for (int i = 0; i < nBin; i ++) {
	 double value = pow(total_true_mass.read(i)-jeans_mass.read(i),2.)/(pow(0.1*total_true_mass.read(i),2)+pow( 0.1*jeans_mass.read(i),2.));
	 chi_sq_test+= value;
	 }
	 chi_sq_test = chi_sq_test/nBin;
	 
	 //printf("New chi_sq = %g, Old chi_sq = %g \n", chi_sq_test, chi_sq[iter-1]);
	 
	 if (chi_sq_test != chi_sq_test){
	 make_new = true;
	 }
	 
	 if (make_new == true) {
	 delta = 1.1*delta;
	 }
	 }
	 
	 //	if (sum_rej> nIter) {
	 //		iter --;
	 //		cin >> testing;
	 //	} //then we take it as a rejection
	 
	 //	else {
	 
	 //cin >> testing;
	 
	 if(chi_sq_test < chi_sq[iter-1]){ // accepts
	 chi_sq[iter] = chi_sq_test;
	 x_0[iter]= x0;
	 x_1[iter]= x1;
	 x_2[iter]= x2;				
	 x_3[iter]= x3;
	 printf(">>>>> Accepting <<<<<<< \n");
	 
	 if (chi_sq_test < chi_sq_ref) {
	 point_a = x0;
	 point_b = x1;
	 point_c = x2;
	 point_d = x3;
	 chi_sq_ref = chi_sq_test;
	 }
	 accept++;
	 }
	 
	 
	 else { //there is a small probability of accepting
	 double r = cat_rng_uniform_double_zero_one();
	 double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
	 printf("r = %g   w = %g \n", r, w);
	 
	 if (r > w) { // i do not accept the new configuration
	 chi_sq[iter] = chi_sq[iter-1];
	 x_0[iter]= x_0[iter-1]; 
	 x_1[iter]= x_1[iter-1]; 
	 x_2[iter]= x_2[iter-1]; 				
	 x_3[iter]= x_3[iter-1]; 
	 
	 rej++;
	 
	 if (rej%nRej == 0) delta = 0.8*delta; //VAMOS TER DE ALTERAR ISTO!!!!!!
	 }
	 
	 if (r<w) {
	 chi_sq[iter] = chi_sq_test;
	 x_0[iter]= x0;
	 x_1[iter]= x1;
	 x_2[iter]= x2;				
	 x_3[iter]= x3;
	 accept++;
	 }
	 
	 //}
	 
	 
	 } // we continue doing the same thing
	 
	 tstop = (double)clock()/CLOCKS_PER_SEC;
	 
	 ttime= tstop-tstart; //ttime is how long your code run
	 
	 printf("It took [%g] time to run one LOOP. N = %d, nIter = %d  nAccept = %d, delta = %g\n", ttime, iter, nIter, accept,delta);
	 
	 iter++;
	 }
	 //*/
	double delta_l = delta;
	double delta_r = delta;
	while (iter < nIter){
		double chi_sq_test;
		printf("New Run\n");
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		
		bool make_new = true;
		int sum_rej = 0;
		while (make_new == true) {
			x0 = x_0[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta_l,delta_l));
			x1 = x_1[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta_l,delta_l));
			//x2 = x_2[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			//x3 = x_3[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			
			//bias = funct_elephant_inside_snake(radius, x0, x1, x2, x3, min_plateau, max_plateau);
			bias = funct_min_to_max(radius, x0, x1, min_plateau, max_plateau);
			total_true_mass = total_mass_from_bias(he_mass, bias);
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (total_true_mass.read(i) < gas_mass.read(i)){
					make_new = true;
					//sum_rej++;
				}
			}
			//if(sum_rej>nIter) break;
			
			dm_density = DensityFromMassSpherical(radius, ArrayMinusArray (total_true_mass, gas_mass));
			
			for (int i = 0; i < nBin-1; i++) {
				if (dm_density.read(i) < dm_density.read(i+1)){
					make_new = true;
					//sum_rej++;
				}
			}
			
			jeans_mass = MassFromJeans(radius, dm_density, ArraySQ(dm_disp_vel), dm_beta);
			
			
			
			chi_sq_test = 0.;
			for (int i = 0; i < nBin; i ++) {
				double value = pow(total_true_mass.read(i)-jeans_mass.read(i),2.)/(pow( 0.1*jeans_mass.read(i),2.));
				chi_sq_test+= value;
			}
			chi_sq_test = chi_sq_test/nBin;
			
			//printf("New chi_sq = %g, Old chi_sq = %g \n", chi_sq_test, chi_sq[iter-1]);
			
			if (chi_sq_test != chi_sq_test){
				make_new = true;
			}
			
			
			if (make_new == true) {
				//sum_rej++;
				//if (sum_rej%30) {
				//	delta_l = 1.1*delta_l;
				//}
				
			}
		}
		
		//	if (sum_rej> nIter) {
		//		iter --;
		//		cin >> testing;
		//	} //then we take it as a rejection
		
		//	else {
		
		//cin >> testing;
		if (chi_sq_test < 1.) {
			printf("bef %g now %g \n", chi_sq[iter-1], chi_sq_test);
			cin >> testing;
		}
		if(chi_sq_test < chi_sq[iter-1]){ // accepts
			chi_sq[iter] = chi_sq_test;
			x_0[iter]= x0;
			x_1[iter]= x1;
			//x_2[iter]= x2;				
			//x_3[iter]= x3;
			printf(">>>>> Accepting <<<<<<< \n");
			
			if (chi_sq_test < chi_sq_ref) {
				point_a = x0;
				point_b = x1;
				//	point_c = x2;
				//	point_d = x3;
				chi_sq_ref = chi_sq_test;
			}
			accept++;
		}
		
		
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
			printf("r = %g   w = %g \n", r, w);
			
			if (r > w) { // i do not accept the new configuration
				chi_sq[iter] = chi_sq[iter-1];
				x_0[iter]= x_0[iter-1]; 
				x_1[iter]= x_1[iter-1]; 
				//x_2[iter]= x_2[iter-1]; 				
				//x_3[iter]= x_3[iter-1]; 
				
				rej++;
				
				if (rej%nRej == 0) delta_l = 0.8*delta_l; //VAMOS TER DE ALTERAR ISTO!!!!!!
			}
			
			if (r<w) {
				chi_sq[iter] = chi_sq_test;
				x_0[iter]= x0;
				x_1[iter]= x1;
				//x_2[iter]= x2;				
				//x_3[iter]= x3;
				accept++;
			}
			
			//}
			
			
		} // we continue doing the same thing
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		printf("It took [%g] time to run one LOOP. LEFT chi_sq = %g N = %d, nIter = %d  nAccept = %d, delta = %g\n", chi_sq[iter], ttime, iter, nIter, accept,delta_l);
		
		iter++;
	}
	for (iter = 0; iter <nIter; iter++) {
		outFile0 << iter << " " << chi_sq[iter] << " " << x_0[iter] << " " << x_1[iter] <<endl;
		
	}
	for (int i = 1; i < nIter; i++) {
		chi_sq[i]= 0.;
	}
	iter = 1;
	while (iter < nIter){
		double chi_sq_test;
		printf("New Run\n");
		//printf("New chi_sq = %g, Old chi_sq = %g \n", chi_sq_test, chi_sq[iter-1]);
		
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		
		bool make_new = true;
		int sum_rej = 0;
		while (make_new == true) {
			//x0 = x_0[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta_l,delta_l));
			//x1 = x_1[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta_l,delta_l));
			x2 = x_2[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta_r,delta_r));
			x3 = x_3[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta_r,delta_r));
			
			//bias = funct_elephant_inside_snake(radius, x0, x1, x2, x3, min_plateau, max_plateau);
			//bias = funct_min_to_max(rkpc, x_0, x_1, min_plateau, max_plateau);
			bias = funct_max_to_min(radius, x2, x3, min_plateau, max_plateau);
			total_true_mass = total_mass_from_bias(he_mass, bias);
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (total_true_mass.read(i) < gas_mass.read(i)){
					make_new = true;
					//sum_rej++;
				}
			}
			//if(sum_rej>nIter) break;
			
			dm_density = DensityFromMassSpherical(radius, ArrayMinusArray (total_true_mass, gas_mass));
			
			for (int i = 0; i < nBin-1; i++) {
				if (dm_density.read(i) < dm_density.read(i+1)){
					make_new = true;
					//sum_rej++;
				}
			}
			
			jeans_mass = MassFromJeans(radius, dm_density, ArraySQ(dm_disp_vel), dm_beta);
			
			
			
			double chi_sq_test = 0.;
			for (int i = 0; i < nBin; i ++) {
				double value = pow(total_true_mass.read(i)-jeans_mass.read(i),2.)/(pow( 0.1*jeans_mass.read(i),2.));
				chi_sq_test+= value;
			}
			chi_sq_test = chi_sq_test/nBin;
			
			//printf("New chi_sq = %g, Old chi_sq = %g \n", chi_sq_test, chi_sq[iter-1]);
			
			if (chi_sq_test != chi_sq_test){
				make_new = true;
			}
			
			//if (make_new == true) {
			//	sum_rej++;
			//	if (sum_rej%30) {
			//		delta_r = 1.1*delta_r;
			//	}
			
			//}
		}
		
		//	if (sum_rej> nIter) {
		//		iter --;
		//		cin >> testing;
		//	} //then we take it as a rejection
		
		//	else {
		
		//cin >> testing;
		
		if(chi_sq_test < chi_sq[iter-1]){ // accepts
			chi_sq[iter] = chi_sq_test;
			//x_0[iter]= x0;
			//x_1[iter]= x1;
			x_2[iter]= x2;				
			x_3[iter]= x3;
			printf(">>>>> Accepting <<<<<<< \n");
			
			if (chi_sq_test < chi_sq_ref) {
				//point_a = x0;
				//point_b = x1;
				point_c = x2;
				point_d = x3;
				chi_sq_ref = chi_sq_test;
			}
			accept++;
		}
		
		
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
			printf("r = %g   w = %g \n", r, w);
			
			if (r > w) { // i do not accept the new configuration
				chi_sq[iter] = chi_sq[iter-1];
				//x_0[iter]= x_0[iter-1]; 
				//x_1[iter]= x_1[iter-1]; 
				x_2[iter]= x_2[iter-1]; 				
				x_3[iter]= x_3[iter-1]; 
				
				rej++;
				
				if (rej%nRej == 0) delta_r = 0.8*delta_r; //VAMOS TER DE ALTERAR ISTO!!!!!!
			}
			
			if (r<w) {
				chi_sq[iter] = chi_sq_test;
				//x_0[iter]= x0;
				//x_1[iter]= x1;
				x_2[iter]= x2;				
				x_3[iter]= x3;
				accept++;
			}
			
			//}
			
			
		} // we continue doing the same thing
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		printf("It took [%g] time to run one LOOP. Right chi_sq = %g N = %d, nIter = %d  nAccept = %d, delta = %g\n", chi_sq[iter], ttime, iter, nIter, accept,delta_r);
		
		iter++;
	}
	
	
	for (iter = 0; iter <nIter; iter++) {
		outFile1 << iter << " " << chi_sq[iter]  <<" " << x_2[iter] << " " << x_3[iter]<<endl;
		
	}
	//outFile.close();
	
	bias.delete_array(); total_true_mass.delete_array(); dm_density.delete_array(); 
	
	
	printf("Best_CHI_SQ = %g for (x0,x1,x2,x3) = (%g,%g,%g,%g)\n", chi_sq_ref, point_a, point_b, point_c, point_d);
	//cin >> testing;
	x_i_return[0] = point_a; x_i_return[1] = point_b; x_i_return[2] = point_c; x_i_return[3] = point_d;
	//cin >> testing;
	
	return x_i_return;
	
}


double * A2052_monte_carlo_loop_step_min_max(int nBin, array radius, array he_mass,  array jeans_mass, 
											 array gas_mass, array dm_beta, array dm_disp_vel, double * x_guess, int nIter, int nRej,double delta, string name_file){
	
	
	array he_mass_cut, gas_mass_cut, jeans_mass_cut, dm_beta_cut, dm_disp_vel_cut, radius_cut;
	he_mass_cut.init(nBin); gas_mass_cut.init(nBin); jeans_mass_cut.init(nBin); dm_beta_cut.init(nBin); dm_disp_vel_cut.init(nBin); radius_cut.init(nBin);
	
	
	for (int i = 0; i < nBin; i++) {
		he_mass_cut.set(i, he_mass.read(i));
		gas_mass_cut.set(i, gas_mass_cut.read(i));
		jeans_mass_cut.set(i, jeans_mass.read(i));
		dm_beta_cut.set(i, dm_beta.read(i));	
		dm_disp_vel_cut.set(i, dm_disp_vel.read(i));	
		radius_cut.set(i, radius.read(i));	
	}
	
	
	array total_mass, dm_density, bias;
	total_mass.init(nBin); dm_density.init(nBin); bias.init(nBin);
	
	
	double * x_return; x_return = new double [2];
	
	double x0_fit[nIter], x1_fit[nIter]; x0_fit[0] = x_guess[0]; x1_fit[0] = x_guess[1];
	double min_plateau = x_guess[2]; double max_plateau = x_guess[3];
	
	double chi_sq[nIter];
	
	//INITIAL CHI_SQ FOR REF
	double chi_sq_test=0.;
	for (int i = 0; i < nBin; i++) {
		double value = pow(he_mass_cut.read(i)-jeans_mass_cut.read(i),2.)/
		(pow(0.1*he_mass_cut.read(i),2.)+ pow(0.1*jeans_mass_cut.read(i),2.));
		chi_sq_test+= value;
	}
	chi_sq_test = chi_sq_test/nBin;
	//printf("chiSQ init = %g \n", chi_sq_test);
	
	chi_sq[0] = chi_sq_test;
	
	// INIT METROPOLIS
	
	double point_a, point_b, chi_ref;
	chi_ref = chi_sq_test;
	
	double tstart, tstop, ttime;
	
	int iter =1; int rej = 0.; int accept = 0;
	double x0_test, x1_test;
	int testing; int sum = 0;
	while (iter<nIter){
		
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		//printf("METROPOLIS: New configuration %d previous chiSQ = %g when (x0,x1) = (%g,%g) delta = %g\n", iter, chi_sq[iter-1], x0_fit[iter-1], x1_fit[iter-1], delta);
		bool make_new = true;
		while (make_new == true) {
			make_new = false;
			x0_test = x0_fit[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			x1_test = x1_fit[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			//printf("LEFT xo, x1 = %g , %g delta = %g \n", x0_test, x1_test, delta);
			//cin >> testing;
			bias = funct_min_to_max(radius_cut, x0_test, x1_test, min_plateau, max_plateau);
			total_mass = total_mass_from_bias(he_mass_cut, bias);
			
			make_new = false;
			for (int j = 0; j < nBin; j++) {
				if (total_mass.read(j) < gas_mass_cut.read(j)) make_new = true;
			}
			
			for (int j = 0; j < nBin-1; j++) {
				if (total_mass.read(j+1)<=total_mass.read(j)) make_new = true;
			}
			dm_density = DensityFromMassSpherical(radius_cut, ArrayMinusArray(total_mass, gas_mass_cut));
			
			for (int j =0 ; j < nBin-1; j++) {
				if (dm_density.read(j)<=dm_density.read(j+1)) make_new = true;
			}
			
			jeans_mass_cut = MassFromJeans(radius_cut, dm_density, ArraySQ(dm_disp_vel_cut), dm_beta_cut);
			
			//	for (int j = 0; j < nBin-1; j++) {
			//		if (jeans_mass_cut.read(j+1)<=jeans_mass_cut.read(j)) make_new = true;
			//		printf("here\n");
			//	}			
			/*	writeArrayToFile(radius_cut, gas_mass_cut, "_gas_mass");
			 writeArrayToFile(radius_cut, total_mass, "_total_mass");
			 writeArrayToFile(radius_cut, jeans_mass_cut, "_jeans_mass");
			 writeArrayToFile(radius_cut, dm_density, "_dm_dens");
			 cin >> testing;
			 //*/	
			
			
			chi_sq_test = 0.;
			for (int j = 0; j < nBin; j++) {
				chi_sq_test+= pow(total_mass.read(j)-jeans_mass_cut.read(j),2.)/
				(pow(0.1*total_mass.read(j),2.)+ pow(0.1*jeans_mass_cut.read(j),2.));
			}
			chi_sq_test = chi_sq_test/nBin;
			
			if (chi_sq_test != chi_sq_test) make_new = true;
			
			if(make_new == true) {
				sum ++;
				if (sum %300) delta = 1.1*delta;
				//		delta = 1.1*delta;
				//		printf("delta = %g\n", delta);
				//	}
			}
			//printf("chi = %g\n", chi_sq_test);
		}
		
		//printf("chi %g chi bef %g \n", chi_sq_test, chi_sq[iter-1]);
		
		if(chi_sq_test < chi_sq[iter-1]){ // accepts
			chi_sq[iter] = chi_sq_test;
			x0_fit[iter]= x0_test;
			x1_fit[iter]= x1_test;
			//printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[iter-1]);
			
			if (chi_sq_test < chi_ref) {
				point_a = x0_test;
				point_b = x1_test;
				chi_ref = chi_sq_test;
			}
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
			//printf("r = %g   w = %g \n", r, w);
			
			if (r > w) { // i do not accept the new configuration
				chi_sq[iter] = chi_sq[iter-1];
				x0_fit[iter]= x0_fit[iter-1]; 
				x1_fit[iter]= x1_fit[iter-1]; 
				
				rej++;
				
				if (rej%nRej == 0) delta = 0.9*delta;
			}
			
			if (r<w) {
				chi_sq[iter] = chi_sq_test;
				x0_fit[iter]= x0_test;
				x1_fit[iter]= x1_test;
				accept++;
			}
			
		}
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		//printf("It took [%g] time to run one LOOP. chi_sq = %g nIter = %d  nAccept = %d\n", ttime, chi_sq_test, nIter, accept);
		
		iter++;
	}
	
	ofstream outFile;
	string MyFile = name_file + "_good_chiSQ_x01.txt";
	outFile.open(MyFile.c_str(), ios::out);
	
	for (int i = 0; i < nIter; i++) {
		outFile << i << " " << chi_sq[i] << " " << x0_fit[i] << " " << x1_fit[i] << endl;
		
	}
	outFile.close();
	
	total_mass.delete_array(); dm_density.delete_array(); dm_disp_vel_cut.delete_array(); dm_beta_cut.delete_array();bias.delete_array();
	jeans_mass_cut.delete_array(); he_mass_cut.delete_array(); gas_mass_cut.delete_array();
	
	printf("minimum chiSq  = %g (x0, x1) = (%g, %g)\n", chi_ref, point_a, point_b);
	//cin >> testing;
	double * value_param; value_param = new double [2];
	value_param[0] = point_a; value_param[1] = point_b;
	
	return value_param;
}

//*/
double * A2052_monte_carlo_loop_step_max_min(int nBin, array radius, array he_mass,  array jeans_mass, 
											 array gas_mass, array dm_beta, array dm_disp_vel, double * x_guess, int nIter, int nRej,double delta, string name_file){
	
	
	array he_mass_cut, gas_mass_cut, jeans_mass_cut, dm_beta_cut, dm_disp_vel_cut, radius_cut;
	he_mass_cut.init(nBin); gas_mass_cut.init(nBin); jeans_mass_cut.init(nBin); dm_beta_cut.init(nBin); dm_disp_vel_cut.init(nBin); radius_cut.init(nBin);
	
	
	for (int i = 0; i < nBin; i++) {
		he_mass_cut.set(i, he_mass.read(i+(8-nBin)));
		gas_mass_cut.set(i, gas_mass_cut.read(i+8-nBin));
		jeans_mass_cut.set(i, jeans_mass.read(i+8-nBin));
		dm_beta_cut.set(i, dm_beta.read(i+8-nBin));	
		dm_disp_vel_cut.set(i, dm_disp_vel.read(i+8-nBin));	
		radius_cut.set(i, radius.read(i+8-nBin));	
	}
	
	
	array total_mass, dm_density, bias;
	total_mass.init(nBin); dm_density.init(nBin); bias.init(nBin);
	
	
	double * x_return; x_return = new double [2];
	
	double x0_fit[nIter], x1_fit[nIter]; x0_fit[0] = x_guess[0]; x1_fit[0] = x_guess[1];
	double min_plateau = x_guess[2]; double max_plateau = x_guess[3];
	
	double chi_sq[nIter];
	
	//INITIAL CHI_SQ FOR REF
	double chi_sq_test=0.;
	for (int i = 0; i < nBin; i++) {
		double value = pow(he_mass_cut.read(i)-jeans_mass_cut.read(i),2.)/(pow(0.1*he_mass_cut.read(i),2.)+ pow(0.1*jeans_mass_cut.read(i),2.));
		chi_sq_test+= value;
	}
	chi_sq_test = chi_sq_test/nBin;
	//printf("chiSQ init = %g \n", chi_sq_test);
	
	chi_sq[0] = chi_sq_test;
	
	// INIT METROPOLIS
	
	double point_a, point_b, chi_ref;
	chi_ref = chi_sq_test;
	
	double tstart, tstop, ttime;
	
	int iter =1; int rej = 0.; int accept = 0;
	double x0_test, x1_test;
	int testing; int sum = 0;
	while (iter<nIter){
		
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		//	printf("METROPOLIS: New configuration %d previous chiSQ = %g when (x0,x1) = (%g,%g) delta = %g\n", iter, chi_sq[iter-1], x0_fit[iter-1], x1_fit[iter-1], delta);
		bool make_new = true;
		while (make_new == true) {
			make_new = false;
			x0_test = x0_fit[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			x1_test = x1_fit[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			//printf("RIGHT xo, x1 = %g , %g delta = %g \n", x0_test, x1_test, delta);
			//cin >> testing;
			bias = funct_max_to_min(radius_cut, x0_test, x1_test, min_plateau, max_plateau);
			total_mass = total_mass_from_bias(he_mass_cut, bias);
			
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (total_mass.read(i) < gas_mass_cut.read(i)) make_new = true;
			}
			
			for (int j = 0; j < nBin-1; j++) {
				if (total_mass.read(j+1)<=total_mass.read(j)) make_new = true;
			}
			dm_density = DensityFromMassSpherical(radius_cut, ArrayMinusArray(total_mass, gas_mass_cut));
			
			for (int j =0 ; j < nBin-1; j++) {
				if (dm_density.read(j)<=dm_density.read(j+1)) make_new = true;
			}
			
			jeans_mass_cut = MassFromJeans(radius_cut, dm_density, ArraySQ(dm_disp_vel_cut), dm_beta_cut);
			
			
			/*	writeArrayToFile(radius_cut, gas_mass_cut, "_gas_mass");
			 writeArrayToFile(radius_cut, total_mass, "_total_mass");
			 writeArrayToFile(radius_cut, jeans_mass_cut, "_jeans_mass");
			 writeArrayToFile(radius_cut, dm_density, "_dm_dens");
			 cin >> testing;
			 //*/	
			
			
			chi_sq_test = 0.;
			for (int j = 0; j < nBin; j++) {
				chi_sq_test+= pow(total_mass.read(j)-jeans_mass_cut.read(j),2.)/(pow(0.1*total_mass.read(j),2.)+ pow(0.1*jeans_mass_cut.read(j),2.));
			}
			chi_sq_test = chi_sq_test/nBin;
			
			if (chi_sq_test != chi_sq_test) make_new = true;
			
			if(make_new == true) {
				sum ++;
				if (sum %300) delta = 1.1*delta;
				//{
				//	delta = 1.1*delta;
				//printf("delta = %g\n", delta);
				//	}
			}
			//printf("chi = %g\n", chi_sq_test);
		}
		
		//printf("chi %g chi bef %g \n", chi_sq_test, chi_sq[iter-1]);
		
		if(chi_sq_test < chi_sq[iter-1]){ // accepts
			chi_sq[iter] = chi_sq_test;
			x0_fit[iter]= x0_test;
			x1_fit[iter]= x1_test;
			//printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[iter-1]);
			
			if (chi_sq_test < chi_ref) {
				point_a = x0_test;
				point_b = x1_test;
				chi_ref = chi_sq_test;
			}
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
			//printf("r = %g   w = %g \n", r, w);
			
			if (r > w) { // i do not accept the new configuration
				chi_sq[iter] = chi_sq[iter-1];
				x0_fit[iter]= x0_fit[iter-1]; 
				x1_fit[iter]= x1_fit[iter-1]; 
				
				rej++;
				
				if (rej%nRej == 0) delta = 0.9*delta;
			}
			
			if (r<w) {
				chi_sq[iter] = chi_sq_test;
				x0_fit[iter]= x0_test;
				x1_fit[iter]= x1_test;
				accept++;
			}
			
		}
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		//	printf("It took [%g] time to run one LOOP iter =%d. chi_sq = %g nIter = %d  nAccept = %d\n", ttime, iter, chi_sq_test, nIter, accept);
		
		iter++;
	}
	
	ofstream outFile;
	string MyFile = name_file + "_good_chiSQ_x23.txt";
	outFile.open(MyFile.c_str(), ios::out);
	
	for (int i = 0; i < nIter; i++) {
		outFile << i << " " << chi_sq[i] << " " << x0_fit[i] << " " << x1_fit[i] << endl;
		
	}
	outFile.close();
	
	
	
	
	total_mass.delete_array(); dm_density.delete_array(); dm_disp_vel_cut.delete_array(); dm_beta_cut.delete_array();bias.delete_array();
	jeans_mass_cut.delete_array(); he_mass_cut.delete_array(); gas_mass_cut.delete_array();
	
	printf("minimum chiSq  = %g (x0, x1) = (%g, %g)\n", chi_ref, point_a, point_b);
	//cin >> testing;
	double * value_param; value_param = new double [2];
	value_param[0] = point_a; value_param[1] = point_b;
	
	return value_param;
}



double * A2052_monte_carlo_loop_step_max_min_left(int nBin, array radius, array he_mass,  array jeans_mass, 
												  array gas_mass, array dm_beta, array dm_disp_vel, double * x_guess, int nIter, int nRej,double delta, string name_file){
	
	
	array he_mass_cut, gas_mass_cut, jeans_mass_cut, dm_beta_cut, dm_disp_vel_cut, radius_cut;
	he_mass_cut.init(nBin); gas_mass_cut.init(nBin); jeans_mass_cut.init(nBin); dm_beta_cut.init(nBin); dm_disp_vel_cut.init(nBin); radius_cut.init(nBin);
	
	
	for (int i = 0; i < nBin; i++) {
		he_mass_cut.set(i, he_mass.read(i));
		gas_mass_cut.set(i, gas_mass_cut.read(i));
		jeans_mass_cut.set(i, jeans_mass.read(i));
		dm_beta_cut.set(i, dm_beta.read(i));	
		dm_disp_vel_cut.set(i, dm_disp_vel.read(i));	
		radius_cut.set(i, radius.read(i));	
	}
	
	array total_mass, dm_density, bias;
	total_mass.init(nBin); dm_density.init(nBin); bias.init(nBin);
	
	
	double * x_return; x_return = new double [2];
	
	double x0_fit[nIter], x1_fit[nIter]; x0_fit[0] = x_guess[0]; x1_fit[0] = x_guess[1];
	double min_plateau = x_guess[2]; double max_plateau = x_guess[3];
	
	double chi_sq[nIter];
	
	//INITIAL CHI_SQ FOR REF
	double chi_sq_test=0.;
	for (int i = 0; i < nBin; i++) {
		double value = pow(he_mass_cut.read(i)-jeans_mass_cut.read(i),2.)/(pow(0.1*he_mass_cut.read(i),2.)+ pow(0.1*jeans_mass_cut.read(i),2.));
		chi_sq_test+= value;
	}
	chi_sq_test = chi_sq_test/nBin;
	//printf("chiSQ init = %g \n", chi_sq_test);
	
	chi_sq[0] = chi_sq_test;
	
	// INIT METROPOLIS
	
	double point_a, point_b, chi_ref;
	chi_ref = chi_sq_test;
	
	double tstart, tstop, ttime;
	
	int iter =1; int rej = 0.; int accept = 0;
	double x0_test, x1_test;
	int testing; int sum = 0;
	while (iter<nIter){
		
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		//	printf("METROPOLIS: New configuration %d previous chiSQ = %g when (x0,x1) = (%g,%g) delta = %g\n", iter, chi_sq[iter-1], x0_fit[iter-1], x1_fit[iter-1], delta);
		bool make_new = true;
		while (make_new == true) {
			make_new = false;
			x0_test = x0_fit[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			x1_test = x1_fit[iter-1] * pow(10.,cat_rng_uniform_double_min_max(-delta,delta));
			//printf("RIGHT xo, x1 = %g , %g delta = %g \n", x0_test, x1_test, delta);
			//cin >> testing;
			bias = funct_max_to_min(radius_cut, x0_test, x1_test, min_plateau, max_plateau);
			total_mass = total_mass_from_bias(he_mass_cut, bias);
			
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (total_mass.read(i) < gas_mass_cut.read(i)) make_new = true;
			}
			
			for (int j = 0; j < nBin-1; j++) {
				if (total_mass.read(j+1)<=total_mass.read(j)) make_new = true;
			}
			dm_density = DensityFromMassSpherical(radius_cut, ArrayMinusArray(total_mass, gas_mass_cut));
			
			for (int j =0 ; j < nBin-1; j++) {
				if (dm_density.read(j)<=dm_density.read(j+1)) make_new = true;
			}
			
			jeans_mass_cut = MassFromJeans(radius_cut, dm_density, ArraySQ(dm_disp_vel_cut), dm_beta_cut);
			
			
			/*	writeArrayToFile(radius_cut, gas_mass_cut, "_gas_mass");
			 writeArrayToFile(radius_cut, total_mass, "_total_mass");
			 writeArrayToFile(radius_cut, jeans_mass_cut, "_jeans_mass");
			 writeArrayToFile(radius_cut, dm_density, "_dm_dens");
			 cin >> testing;
			 //*/	
			
			
			chi_sq_test = 0.;
			for (int j = 0; j < nBin; j++) {
				chi_sq_test+= pow(total_mass.read(j)-jeans_mass_cut.read(j),2.)/(pow(0.1*total_mass.read(j),2.)+ pow(0.1*jeans_mass_cut.read(j),2.));
			}
			chi_sq_test = chi_sq_test/nBin;
			
			if (chi_sq_test != chi_sq_test) make_new = true;
			
			if(make_new == true) {
				sum ++;
				if (sum %300) delta = 1.1*delta;
				//{
				//	delta = 1.1*delta;
				//printf("delta = %g\n", delta);
				//	}
			}
			//printf("chi = %g\n", chi_sq_test);
		}
		
		//printf("chi %g chi bef %g \n", chi_sq_test, chi_sq[iter-1]);
		
		if(chi_sq_test < chi_sq[iter-1]){ // accepts
			chi_sq[iter] = chi_sq_test;
			x0_fit[iter]= x0_test;
			x1_fit[iter]= x1_test;
			//printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[iter-1]);
			
			if (chi_sq_test < chi_ref) {
				point_a = x0_test;
				point_b = x1_test;
				chi_ref = chi_sq_test;
			}
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
			//printf("r = %g   w = %g \n", r, w);
			
			if (r > w) { // i do not accept the new configuration
				chi_sq[iter] = chi_sq[iter-1];
				x0_fit[iter]= x0_fit[iter-1]; 
				x1_fit[iter]= x1_fit[iter-1]; 
				
				rej++;
				
				if (rej%nRej == 0) delta = 0.9*delta;
			}
			
			if (r<w) {
				chi_sq[iter] = chi_sq_test;
				x0_fit[iter]= x0_test;
				x1_fit[iter]= x1_test;
				accept++;
			}
			
		}
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		//	printf("It took [%g] time to run one LOOP iter =%d. chi_sq = %g nIter = %d  nAccept = %d\n", ttime, iter, chi_sq_test, nIter, accept);
		
		iter++;
	}
	
	ofstream outFile;
	string MyFile = name_file + "_good_chiSQ_x23.txt";
	outFile.open(MyFile.c_str(), ios::out);
	
	for (int i = 0; i < nIter; i++) {
		outFile << i << " " << chi_sq[i] << " " << x0_fit[i] << " " << x1_fit[i] << endl;
		
	}
	outFile.close();
	
	
	
	
	total_mass.delete_array(); dm_density.delete_array(); dm_disp_vel_cut.delete_array(); dm_beta_cut.delete_array();bias.delete_array();
	jeans_mass_cut.delete_array(); he_mass_cut.delete_array(); gas_mass_cut.delete_array();
	
	printf("minimum chiSq  = %g (x0, x1) = (%g, %g)\n", chi_ref, point_a, point_b);
	//cin >> testing;
	double * value_param; value_param = new double [2];
	value_param[0] = point_a; value_param[1] = point_b;
	
	return value_param;
}



array metropolis_given_all(array radius, array  he_mass,  array jeans_mass, array gas_mass, array dm_beta, array dm_disp_vel, array bias_ver, int nIter, int  nRej, string name_file){
	
	int nBin = radius.numberBins();
	
	array total_mass, dm_density, bias;
	total_mass.init(nBin); dm_density.init(nBin); bias.init(nBin);
	
	
	double bias_test; double bias_fit[nIter];
	bias_fit[0] = 0.;
	
	double delta = .99;
	double chi_sq[nIter];
	
	//INITIAL CHI_SQ FOR REF
	double chi_sq_test=0.;
	for (int i = 0; i < nBin; i++) {
		double value = pow(he_mass.read(i)-jeans_mass.read(i),2.)/(pow(0.1*he_mass.read(i),2.)+ pow(0.1*jeans_mass.read(i),2.));
		chi_sq_test+= value;
	}
	chi_sq_test = chi_sq_test/nBin;
	//printf("chiSQ init = %g \n", chi_sq_test);
	
	chi_sq[0] = chi_sq_test;
	
	// INIT METROPOLIS
	
	double point_a, point_b, chi_ref;
	chi_ref = chi_sq_test;
	
	double tstart, tstop, ttime;
	
	int iter =1; int rej = 0.; int accept = 0;
	double x0_test, x1_test;
	int testing; int sum = 0;
	while (iter<nIter){
		
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		//	printf("METROPOLIS: New configuration %d previous chiSQ = %g when (x0,x1) = (%g,%g) delta = %g\n", iter, chi_sq[iter-1], x0_fit[iter-1], x1_fit[iter-1], delta);
		bool make_new = true;
		while (make_new == true) {
			make_new = false;
			bias_test = bias_fit[iter-1] + cat_rng_uniform_double_min_max(-delta,delta);
			
			for (int i = 0; i < nBin; i++) {
				if(i==0)bias.set(0, bias_test);
				else bias.set(i, bias_ver.read(i));
			}
			
			
			
			//bias = funct_max_to_min(radius_cut, x0_test, x1_test, min_plateau, max_plateau);
			total_mass = total_mass_from_bias(he_mass, bias);
			
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (total_mass.read(i) < gas_mass.read(i)) make_new = true;
			}
			
			for (int j = 0; j < nBin-1; j++) {
				if (total_mass.read(j+1)<=total_mass.read(j)) make_new = true;
			}
			dm_density = DensityFromMassSpherical(radius, ArrayMinusArray(total_mass, gas_mass));
			
			for (int j =0 ; j < nBin-1; j++) {
				if (dm_density.read(j)<=dm_density.read(j+1)) make_new = true;
			}
			
			jeans_mass = MassFromJeans(radius, dm_density, ArraySQ(dm_disp_vel), dm_beta);
			
			
			/*	writeArrayToFile(radius_cut, gas_mass_cut, "_gas_mass");
			 writeArrayToFile(radius_cut, total_mass, "_total_mass");
			 writeArrayToFile(radius_cut, jeans_mass_cut, "_jeans_mass");
			 writeArrayToFile(radius_cut, dm_density, "_dm_dens");
			 cin >> testing;
			 //*/	
			
			
			chi_sq_test = 0.;
			for (int j = 0; j < nBin; j++) {
				chi_sq_test+= pow(total_mass.read(j)-jeans_mass.read(j),2.)/(pow(0.1*total_mass.read(j),2.)+ pow(0.1*jeans_mass.read(j),2.));
			}
			chi_sq_test = chi_sq_test/nBin;
			
			if (chi_sq_test != chi_sq_test) make_new = true;
			
			if(make_new == true) {
				sum ++;
				if (sum %300) delta = 1.1*delta;
				//{
				//	delta = 1.1*delta;
				//printf("delta = %g\n", delta);
				//	}
			}
			//printf("chi = %g\n", chi_sq_test);
		}
		
		//printf("chi %g chi bef %g \n", chi_sq_test, chi_sq[iter-1]);
		
		if(chi_sq_test < chi_sq[iter-1]){ // accepts
			chi_sq[iter] = chi_sq_test;
			bias_fit[iter]= bias_test;
			//x1_fit[iter]= x1_test;
			//printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[iter-1]);
			
			if (chi_sq_test < chi_ref) {
				point_a = bias_test;
				//point_b = x1_test;
				chi_ref = chi_sq_test;
			}
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
			//printf("r = %g   w = %g \n", r, w);
			
			if (r > w) { // i do not accept the new configuration
				chi_sq[iter] = chi_sq[iter-1];
				bias_fit[iter]= bias_fit[iter-1]; 
				//x1_fit[iter]= x1_fit[iter-1]; 
				
				rej++;
				
				if (rej%nRej == 0) delta = 0.9*delta;
			}
			
			if (r<w) {
				chi_sq[iter] = chi_sq_test;
				bias_fit[iter]= bias_test;
				//x1_fit[iter]= x1_test;
				accept++;
			}
			
		}
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		//	printf("It took [%g] time to run one LOOP iter =%d. chi_sq = %g nIter = %d  nAccept = %d\n", ttime, iter, chi_sq_test, nIter, accept);
		
		iter++;
	}
	
	ofstream outFile;
	string MyFile = name_file + "_good_chiSQ_bias.txt";
	outFile.open(MyFile.c_str(), ios::out);
	
	for (int i = 0; i < nIter; i++) {
		outFile << i << " " << chi_sq[i] << " " << bias_fit[i]  << endl;
		
	}
	outFile.close();
	
	
	
	
	total_mass.delete_array(); dm_density.delete_array();
	
	printf("minimum chiSq  = %g (x0, x1) = (%g, %g)\n", chi_ref, point_a, point_b);
	//cin >> testing;
	
	for (int i = 0; i < nBin; i++) 
	{
		if(i==0)bias.set(0, point_a);
		else bias.set(i, bias_ver.read(i));
	}
	
	return bias;
}


cluster init_cluster_with_limits( double R200_kpc, int nFit, double r_min_R200, double r_max_R200){
	
	int n_fit = nFit;
	
	cluster cluster_fitted;
	cluster_fitted.SetNConfig(100);//it does not matter just to begin
	cluster_fitted.SetNBin(n_fit);
	
	cluster_fitted.SetR200(R200_kpc);
	cluster_fitted.SetRMin(r_min_R200*R200_kpc*1.e3*onePcMtr, 1.e17);	// in meters
	cluster_fitted.SetRMax(r_max_R200*R200_kpc*1.e3*onePcMtr, 1.e22);	// in meters
	cluster_fitted.InitCluster();
	
	//for(int i = 0; i < n_fit; i++) printf("r=%g, r200+=%g\n", cluster_fitted.RLin[0].read(i), cluster_fitted.Rkpc[0].read(i));
	
	return cluster_fitted;	
}
//initializes another array this time with narrow radial band
double * fit_double_power_law_return_fit_and_error(array x_data, array y_data, double * guess_parameters){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = .1*y_fitting[i];
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double [8];
	fit = cat_fit_doublePL(guess_parameters, &my_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return fit;
	
}

double * fit_double_power_law_return_fit(array x_data, array y_data, double * guess_parameters){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = 1.*y_data.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[8];
	fit = cat_fit_doublePL(guess_parameters, &my_data); //routine to find the best-fiting parameters that describe data as following a power law

	return fit;
}

double * fit_double_power_law_with_error_return_fit_and_error(array x_data, array y_data, array y_error, double * guess_parameters){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[8];
	fit = cat_fit_doublePL(guess_parameters, &my_data); //routine to find the best-fiting parameters that describe data as following a power law

	//array * data_fitted; data_fitted = new array[2];
	//data_fitted[0].init(nBin); data_fitted[1].init(nBin);
	//data_fitted[0] = DensityPowerLaw(x_data, fit[0], fit[1], fit[2], fit[3]);
	return fit;
	
}
array * fit_obs_to_doublePowerLaw(array cluster_obs_radius, array cluster_obs_data, array cluster_fit_radius, double * guess){
	
	int n_fit = cluster_fit_radius.numberBins();
	
	double * fit; fit = new double[8];
	fit = fit_double_power_law_return_fit_and_error(cluster_obs_radius, cluster_obs_data, guess);
	
	array * doublePL; doublePL = new array [2];
	doublePL[0].init(n_fit); doublePL[1].init(n_fit);
	
	doublePL[0] = DensityPowerLaw(cluster_fit_radius, fit[0], fit[1], fit[2], fit[3]);
	//rho0=fit[0], r_s=fit[1], alpha=fit[2], beta=fit[3]
	
	//	gas_density[1] = ERRO; - falta fazer
	doublePL[1] = DoubleTimesArray(0.1,doublePL[0]);
	
	delete [] fit;
	
	return doublePL;	
}

double * fit_double_power_law_fix_beta(array x_data, array y_data, array y_error, double * guess_parameters){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[8];
	fit = cat_fit_doublePL_fix_beta(guess_parameters, &my_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return fit;
	
}

double * fit_double_power_law_fix_alpha_beta(array x_data, array y_data, array y_error, double * guess_parameters){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[8];
	fit = cat_fit_doublePL_fix_alpha_beta(guess_parameters, &my_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return fit;
	
}

double * fit_parametrisation(array x_data, array y_data, array y_error, double * guess_parameters){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[4];
	fit = cat_fit_parametrisation_temp(guess_parameters, &my_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return fit;
	
}

array * fit_mass_from_density_PL(array dataX, array dataY, double * guess){
	int nBin = dataY.numberBins();
	array * FIT; FIT = new array [2];
	FIT[0].init(nBin);
	FIT[1].init(8);
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = .1*y_dm_mass[i];
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	double rho0, r_s, alpha, beta, err_rho0, err_r_s, err_alpha, err_beta;
	rho0 = p_fit_dm_mass[0]; r_s = p_fit_dm_mass[1]; alpha = p_fit_dm_mass[2]; beta = p_fit_dm_mass[3];
	err_rho0 = p_fit_dm_mass[4];  err_r_s = p_fit_dm_mass[5];  err_alpha = p_fit_dm_mass[6];  err_beta = p_fit_dm_mass[7];
	
	for (int i = 0; i < nBin; i++) y_dm_mass[i] = 4.*pi*rho0*pow(r[i],2.) * pow((r[i]+r_s)/r[i], alpha) * pow(r_s/(r[i] + r_s), beta);
	
	double *dm;
	dm = new double[nBin];
	dm = integralTrapezeZeroToPointP(nBin, r,y_dm_mass);
	
	for (int i = 0; i < nBin; i++) FIT[0].set(i, dm[i]);
	
	for (int i = 0; i < 8; i++) FIT[1].set(i, p_fit_dm_mass[i]);
	//err_mass_assume_density_powerLaw(p_fit_dm_mass, dataX);
	delete[] p_fit_dm_mass;
	delete[] dm;
	
	printf("Sucess - Fit Dm Mass \n");
	
	return FIT;
}
double * fit_mass_from_density_PL_return_parameters(array dataX, array dataY, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = .1*y_dm_mass[i];
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
}
double * fit_mass_from_density_PL_return_parameters_set_alpha(array dataX, array dataY, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = .1*y_dm_mass[i];
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M_fix_alpha(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
}


double * fit_mass_from_density_PL_return_parameters_set_r_s_beta(array dataX, array dataY, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = .1*y_dm_mass[i];
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M_fix_r_s_beta(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
}
double * fit_mass_from_density_PL_return_parameters_w_err(array dataX, array dataY, array err, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = err.read(i);
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
}

														  
double * fit_mass_from_density_PL_return_parameters_set_alpha_beta(array dataX, array dataY, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = .1*y_dm_mass[i];
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M_fix_alpha_beta(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
}

double * fit_density_PL_w_error_return_parameters_set_alpha_beta(array dataX, array dataY, array dataY_error, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = dataY_error.read(i);
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_fix_alpha_beta(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
}

double * fit_density_PL_w_error_return_parameters_set_rs(array dataX, array dataY, array dataY_error, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = dataY_error.read(i);
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_fix_rs(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
}


double * fit_mass_from_density_Sersic_return_parameters(array dataX, array dataY, double * guess){
	int nBin = dataY.numberBins();
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = dataX.read(i);
		y_dm_mass[i] = dataY.read(i);
		sigma_dm_mass[i] = .1*y_dm_mass[i];
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M_fix_alpha_beta(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	return p_fit_dm_mass;
} //not done
														






double * fit_beta_profile_with_error_return_fit_and_error(array x_data, array y_data, array y_error, double * beta_guess){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[6];
	fit = cat_fit_beta_profile(beta_guess,&my_data); //routine to find the best-fiting parameters that describe data as following a power law
	//array * data_fitted; data_fitted = new array[2];
	//data_fitted[0].init(nBin); data_fitted[1].init(nBin);
	//data_fitted[0] = DensityPowerLaw(x_data, fit[0], fit[1], fit[2], fit[3]);
	//array * data_fitted; data_fitted = new array[1];
	//data_fitted[0].init(nBin);
	//data_fitted[0] = DensityBetaModel(x_data, fit[0], fit[1], fit[2]);
	//data_fitted[1] = aqui poriamos o erro
	//delete [] fit;
	return fit;
	
}
double * fit_sersic_profile_with_error_return_fit_and_error(array x_data, array y_data, array y_error, double * sersic_guess){
	printf("Making Sersic model\n");
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[6];
	fit = cat_fit_sersic_profile(sersic_guess,&my_data); //routine to find the best-fiting parameters that describe data as following a power law
	//array * data_fitted; data_fitted = new array[2];
	//data_fitted[0].init(nBin); data_fitted[1].init(nBin);
	//data_fitted[0] = DensityPowerLaw(x_data, fit[0], fit[1], fit[2], fit[3]);
	//array * data_fitted; data_fitted = new array[1];
	//data_fitted[0].init(nBin);
	//data_fitted[0] = DensityBetaModel(x_data, fit[0], fit[1], fit[2]);
	//data_fitted[1] = aqui poriamos o erro
	//delete [] fit;
	return fit;
	
}

double * fit_quadratic_with_error_return_fit_and_error(array x_data, array y_data, array y_error, double * quadratic_guess){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[6];
	fit = cat_fit_quadratic(quadratic_guess, &my_data); //routine to find the best-fiting parameters that describe data as following a power law
	//array * data_fitted; data_fitted = new array[2];
	//data_fitted[0].init(nBin); data_fitted[1].init(nBin);
	//data_fitted[0] = DensityPowerLaw(x_data, fit[0], fit[1], fit[2], fit[3]);
	//array * data_fitted; data_fitted = new array[1];
	//data_fitted[0].init(nBin);
	//data_fitted[0] = DensityBetaModel(x_data, fit[0], fit[1], fit[2]);
	//data_fitted[1] = aqui poriamos o erro
	//delete [] fit;
	return fit;
	
}
array * fit_obs_to_quadratic(array cluster_obs_radius, array cluster_obs_data, array cluster_fit_radius, double * guess){
	
	int n_fit = cluster_fit_radius.numberBins();
	
	double * fit; fit = new double[6];
	fit = fit_quadratic_with_error_return_fit_and_error(cluster_obs_radius, cluster_obs_data, DoubleTimesArray(0.1,cluster_obs_data),guess);
	
	array * quadratic; quadratic = new array [2];
	quadratic[0].init(n_fit); quadratic[1].init(n_fit);
	
	for (int i = 0; i < n_fit; i++) {
		quadratic[0].set(i,fit[0]*pow(cluster_fit_radius.read(i),2.) + fit[1]*cluster_fit_radius.read(i) + fit[2]);
	}
	
	//doublePL[0] = DensityPowerLaw(cluster_fit_radius, fit[0], fit[1], fit[2], fit[3]);
	//rho0=fit[0], r_s=fit[1], alpha=fit[2], beta=fit[3]
	
	//	gas_density[1] = ERRO; - falta fazer
	quadratic[1] = DoubleTimesArray(0.1,quadratic[0]);
	
	return quadratic;	
}

double * fit_exp_decay_with_error_return_fit_and_error(array x_data, array y_data, array y_error, double * quadratic_guess){
	
	int nBin = x_data.numberBins();
	
	double r[nBin], y_fitting[nBin], sigma_y[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = x_data.read(i);
		y_fitting[i] = y_data.read(i);
		sigma_y[i] = y_error.read(i);
	}
	struct mini_cluster my_data = {nBin, r, y_fitting, sigma_y};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * fit; fit = new double[6];
	fit = cat_fit_exp(3, nBin, sigma_y, y_fitting);
	//fit = cat_fit_quadratic(quadratic_guess,3, &my_data); //routine to find the best-fiting parameters that describe data as following a power law
	//array * data_fitted; data_fitted = new array[2];
	//data_fitted[0].init(nBin); data_fitted[1].init(nBin);
	//data_fitted[0] = DensityPowerLaw(x_data, fit[0], fit[1], fit[2], fit[3]);
	//array * data_fitted; data_fitted = new array[1];
	//data_fitted[0].init(nBin);
	//data_fitted[0] = DensityBetaModel(x_data, fit[0], fit[1], fit[2]);
	//data_fitted[1] = aqui poriamos o erro
	//delete [] fit;
	return fit;
	
}
array * fit_obs_to_exp_decay(array cluster_obs_radius, array cluster_obs_data, array cluster_fit_radius, double * guess){
	
	int n_fit = cluster_fit_radius.numberBins();
	
	double * fit; fit = new double[6];
	fit = fit_exp_decay_with_error_return_fit_and_error(cluster_obs_radius, cluster_obs_data, DoubleTimesArray(0.1,cluster_obs_data),guess);
	
	array * quadratic; quadratic = new array [2];
	quadratic[0].init(n_fit); quadratic[1].init(n_fit);
	
	for (int i = 0; i < n_fit; i++) {
		double value = fit[0] * exp(-fit[1] * i) + fit[2] ;
		quadratic[0].set(i,value);
	}
	
	//doublePL[0] = DensityPowerLaw(cluster_fit_radius, fit[0], fit[1], fit[2], fit[3]);
	//rho0=fit[0], r_s=fit[1], alpha=fit[2], beta=fit[3]
	
	//	gas_density[1] = ERRO; - falta fazer
	quadratic[1] = DoubleTimesArray(0.1,quadratic[0]);
	
	return quadratic;	
}

array funct_min_to_max(array r_R200, double r0_R200, double r1_R200, double min_plateau, double max_plateau){

	int nBin = r_R200.numberBins();
	
	array function; function.init(nBin);
	
//	yd(x) = mlog(xd)+b
//	ye(x) = mlog(xe)+b
//	yd + ye = mlog(xd*xe) + 2b
//	yd - ye = mlog(xd/xe)
	
	double m = (max_plateau - min_plateau)/ log10(r1_R200/r0_R200);
	double b = .5 * (max_plateau + min_plateau - m*log10(r1_R200*r0_R200));
					 
	for (int i = 0; i<nBin; i++) {
		double r = r_R200.read(i);
		double value = 0.;
		if(r < r0_R200) value = min_plateau;
		
		if(r > r1_R200) value = max_plateau;
		
		if (r >= r0_R200 && r <= r1_R200) value = m * log10(r) + b;
	
		function.set(i, value);
	}
	
	return function;
					 
}
array funct_max_to_min(array r_R200, double r0_R200, double r1_R200, double min_plateau, double max_plateau){
	
	int nBin = r_R200.numberBins();
	
	array function; function.init(nBin);
	
	//	yd(x) = mlog(xd)+b
	//	ye(x) = mlog(xe)+b
	//	yd + ye = mlog(xd*xe) + 2b
	//	yd - ye = mlog(xd/xe)
	
	double m = -(max_plateau - min_plateau)/ log10(r1_R200/r0_R200);
	double b = .5 * (max_plateau + min_plateau - m*log10(r1_R200*r0_R200));
	double value = 0.;
	
	for (int i = 0; i<nBin; i++) {
		double r = r_R200.read(i);

		if(r < r0_R200) value = max_plateau;
		
		if(r > r1_R200) value = min_plateau;
		
		if (r >= r0_R200 && r <= r1_R200) value = m * log10(r) + b;
		//printf("i= %d and value = %g \n", i, value);
		function.set(i, value);
	}
	
	return function;
	
}
array funct_elephant_inside_snake(array r_R200, double r0l_R200, double r1l_R200,double r0r_R200, double r1r_R200, double min_plateau, double max_plateau){
	
	int nBin = r_R200.numberBins();
	
	array function; function.init(nBin);
	
	function = ArrayPlusArray(funct_min_to_max(r_R200, r0l_R200, r1l_R200, min_plateau, max_plateau),
							  funct_max_to_min(r_R200, r0r_R200, r1r_R200, min_plateau, max_plateau));
	
	return function;
	
}

array new_funct_biased_elephant_inside_snake (array radius_obs, double x_left_0, double x_left_1, int n_left, 
											  double x_right_0, double x_right_1, int n_right, double min_plateau, double max_plateau){
	
	int nBin = radius_obs.numberBins();
	
	array function; function.init(nBin);
	array f_l, f_r; f_l.init(nBin); f_r.init(nBin);
	f_l = funct_min_to_max(radius_obs, x_left_0,x_left_1,  min_plateau,  max_plateau);
	f_r = funct_max_to_min(radius_obs, x_right_0,x_right_1,  min_plateau,  max_plateau);
	
	for (int i = 0; i < nBin;i++) {
		if (i < n_left) function.set(i, f_l.read(i));
		else function.set(i, 0.);
	}
	for (int i = 0; i < nBin;i++) {
		if (i > (nBin-n_right)) function.set(i, function.read(i)+ f_r.read(i));
	}
	
	f_l.delete_array(); f_r.delete_array();
	return function;
	
}

array total_mass_from_bias_elephant_inside_snake(array mass_he, array r_R200, double r0l_R200, double r1l_R200,double r0r_R200, double r1r_R200, double min_plateau, double max_plateau){
	
	int nBin = r_R200.numberBins();
	
	array total_mass; total_mass.init(nBin);
	
	total_mass = ArrayDividesArray(mass_he, 
								   ArrayPlusInt(funct_elephant_inside_snake(r_R200,r0l_R200,r1l_R200,
																			r0r_R200,r1r_R200,min_plateau,max_plateau),
												1));

	return total_mass;
	
}
array total_mass_from_bias(array mass_he, array bias){
	int nBin = mass_he.numberBins();
	
	array total_mass; total_mass.init(nBin);
	
	total_mass = ArrayDividesArray(mass_he, 
								   ArrayPlusInt(bias,1));
	
	return total_mass;
	
}





void cluster::Fit_dm_mass(double * guess){
	//saves Cluster.dm_mass[0]
	//saves Cluster.dm_mass_error[0]
	
	
	int nBin = dm_mass[0].numberBins();
	double rho0, r_s, alpha, beta;
	double err_rho0, err_r_s, err_alpha, err_beta;
	
	double r[nBin], y_dm_mass[nBin], sigma_dm_mass[nBin]; 
	for (int i = 0; i < nBin; i++) {
		r[i] = RLin[0].read(i);
		y_dm_mass[i] = dm_mass[0].read(i);
		sigma_dm_mass[i] = .1*y_dm_mass[i];
	}
	struct mini_cluster dm_data = {nBin, r, y_dm_mass, sigma_dm_mass};  //assignes a structure for dm_data: radius, mass and std deviation
	
	double * p_fit_dm_mass; p_fit_dm_mass[8];
	p_fit_dm_mass = cat_fit_doublePL_M(guess,&dm_data); //routine to find the best-fiting parameters that describe data as following a power law
	
	rho0 = p_fit_dm_mass[0]; r_s = p_fit_dm_mass[1]; alpha = p_fit_dm_mass[2]; beta = p_fit_dm_mass[3];
	err_rho0 = p_fit_dm_mass[4];  err_r_s = p_fit_dm_mass[5];  err_alpha = p_fit_dm_mass[6];  err_beta = p_fit_dm_mass[7];
	
	for (int i = 0; i < 8; i++) parameters_fit_dm_mass[0].set(i, p_fit_dm_mass[i]);
	
	for (int i = 0; i < nBin; i++) y_dm_mass[i] = 4.*pi*rho0*pow(r[i],2.) * pow((r[i]+r_s)/r[i], alpha) * pow(r_s/(r[i] + r_s), beta);
	
	double *dm;
	dm = new double[nBin];
	dm = integralTrapezeZeroToPointP(nBin, r,y_dm_mass);
	
	for (int i = 0; i < nBin; i++) dm_mass[0].set(i, dm[i]);

	Err_dm_mass_assume_density_powerLaw(nBin); //saves dm_mass_error = ∆(M_dm)
	
	printf("Sucess - Fit Dm Mass \n");
	delete[] p_fit_dm_mass;
	delete[] dm;
	
}
void cluster::Init_dm_properties_beta_zero(double * guess){
	
	int nBin = RLin[0].numberBins();
	
	mass_HE[0] = MassFromHE(RLin[0],gas_density[0],gas_temperature[0]);
	printf("Calculated Mass He from Gas density and temperature \n");
	
	dm_mass[0] = ArrayMinusArray(mass_HE[0], gas_mass[0]);
	printf("First assuming that Mass He = M total - calculates dm Mass \n");
	
	array * fitting; fitting = new array[2];
	fitting = fit_mass_from_density_PL(RLin[0], dm_mass[0], guess); //return fit for dm_mass +/- ∆ dm_mass
	printf("Since dm Mass will be scattered we fitted this values for a dm Density profile:\n");
	printf("rho0 = %g +/- %g; r_s = %g +/- %g; alpha = %g +/- %g; beta = %g +/- %g\n", 
		   fitting[1].read(0), fitting[1].read(4), fitting[1].read(1), fitting[1].read(5),
		   fitting[1].read(2), fitting[1].read(6), fitting[1].read(3), fitting[1].read(7));
	
	for(int i = 0; i < nBin; i++) dm_mass[0].set(i, fitting[0].read(i)); 
	//dm_density[0] = DensityFromMassSpherical(RLin[0], dm_mass[0]); //from real data  
	dm_density[0] = DensityPowerLaw(RLin[0], 1.e-23*10., 500.*1.e3*onePcMtr, 1., 3.); //from fake data

	//Err_dm_mass_assume_density_powerLaw(nBin);
	
	dm_density_error[0] = err_one_variable (dm_density[0], dm_mass[0],dm_mass_error[0]); //saves dm_density_error = ∆(dens_dm)

	dm_gamma[0] = Derivative(Log10Array(RLin[0]), Log10Array(dm_density[0])); //gamma
	dm_gamma_error[0] = err_one_variable(dm_gamma[0],dm_density[0], dm_density_error[0]);
	
//	dm_beta[0] = BetaFollowsGammaLinear(dm_gamma[0]); //beta from real
	
	dm_beta[0] =  SetAllToZero(nBin); //I took data when beta = 0 let's leave it like this for now
//	Err_dm_beta_from_gamma();
	
	dm_disp_velocity[0] = dispVelFromTemperature(gas_temperature[0], dm_beta[0]); //disp_velocity
	dm_disp_velocity_error[0] = err_one_variable(dm_disp_velocity[0], dm_beta[0], dm_beta_error[0]); //because we still do not have ∆T
	
	mass_jeans[0] = MassFromJeans(RLin[0], dm_density[0], ArraySQ(dm_disp_velocity[0]), dm_beta[0]);
		
	writeArrayToFile(Rkpc[0], ArrayAbs(ArrayDividesArray(ArrayMinusArray(mass_HE[0], 
																		 mass_jeans[0]), mass_HE[0])), "A_fake");
	writeArrayToFile(Rkpc[0], ArrayAbs(ArrayDividesArray( ArrayMinusArray(mass_HE[0], mass_jeans[0])
														 ,ArrayPlusArray(mass_jeans[0], mass_HE[0]))), "A+B_fake");
	
	
	printf("Success - Init Dm Properties \n");
	
	delete[] fitting;
}
void cluster::Non_thermal_mass_init(double r_s){
	
	//REVER PARA FAZER OS GRAFIICOS DA MASSA
	int nBin = gas_density[0].numberBins();
	
	double cr0, b0, exp_cr, exp_b;
	cr0 = 0.46;
	exp_cr = -0.7;
	b0 = 4.6*1.e-6*Gauss;
	exp_b = 0.6;
	
	array cosmic_rays, turbulence, magnetic;
	cosmic_rays.init(nBin);
	turbulence.init(nBin);
	magnetic.init(nBin);
	
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(mass_HE[0],oneMsunkg), "heMass");
	
	SetCRsPressure(cr0, exp_cr);
	YpCRs[0] = DoubleTimesArray(cr0,ArrayPowerDouble(ArrayDividesDouble(RLin[0], r_s), exp_cr));
	cosmic_rays = MassFromCRsPressure(RLin[0], gas_temperature[0], YpCRs[0], mass_HE[0]);
	
	//writeArrayToFile(Rkpc[0], ArrayDividesDouble(cosmic_rays,oneMsunkg), "cosmicMass");
	//for (int i = 0; i < nBin; i++) {
	//	printf("Bias = %g for cosmic_mass = %g and mass_HE = %g\n", 100*(cosmic_rays.read(i))/(cosmic_rays.read(i)+mass_HE[0].read(i)),
	//																cosmic_rays.read(i), mass_HE[0].read(i));
	//}
	//writeArrayToFile(Rkpc[0], ArrayDividesArray(DoubleTimesArray(100.,cosmic_rays), ArrayPlusArray(cosmic_rays, mass_HE[0])), "bias_CRs");
	//SetTurbPressure("case1",false); //This nee
	//TotalMassCalculus("Turb");
	//turbulence = MassFromTurbPressure(RLin[0], v_r[0], v_theta[0]);
	
	SetMagneticParameters(b0, exp_b); //minimum bfield in GC
	int bin;
	int i =0;
	while(i < nBin) {
		//comparison with R200 to get the bin R200 = 1.8e6pc = 1.8Mpc
		if (RLin[0].read(i) > r_s) {
			bin = i;
			//bin = 0.1*nBin;
			printf("Found the BIN %d\n", bin);
			break;
		}
		i++;
	}
	MagneticField[0] = DoubleTimesArray(b0,ArrayPowerDouble(ArrayDividesDouble(gas_density[0],
																			gas_density[0].read(bin)), exp_b));
	magnetic = MassFromMagPressure(RLin[0], gas_density[0], MagneticField[0]); 
	//for (int i = 0; i < nBin; i++) {
	//	cout << i << endl;
	//	printf("Bias = %g for magnetic_mass = %g and mass_HE = %g\n", 100*(-magnetic.read(i))/(magnetic.read(i)+mass_HE[0].read(i)),
		   															//magnetic.read(i), mass_HE[0].read(i));
	//}
		   
	writeArrayToFile(Rkpc[0], ArrayDividesArray(DoubleTimesArray(-100.,magnetic), ArrayPlusArray(magnetic, mass_HE[0])), "magneticMass");
	
	
	//writeArrayToFile(Cluster.Rkpc[0], ArrayDividesDouble(magnetic,oneMsunkg), "massHE");
	//nonThermalMass = ArrayPlusArray( magnetic, ArrayPlusArray(turbulence, cosmic_rays));
	mass_nTh[0] = ArrayPlusArray( magnetic, cosmic_rays);
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(mass_nTh[0],oneMsunkg), "nonThermalMass");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(cosmic_rays,oneMsunkg), "cosmicMass");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(magnetic,oneMsunkg), "masgneticMass");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(ArrayPlusArray(mass_nTh[0], mass_HE[0]),oneMsunkg), "totalMAss");
	
	cosmic_rays.delete_array(); turbulence.delete_array(); magnetic.delete_array();
	//printf("Sucess - First Non thermal mass \n");
}
void cluster::Jeans_mass_from_total_mass_true(array TotalMass, array GasMass_Fit, double * guess){
	
	dm_mass[0] = ArrayMinusArray(TotalMass, GasMass_Fit);
	
	//fit_dm_mass(guess); //saves Cluster.dm_mass[0]
								//saves Cluster.dm_mass_error[0]
	
	printf("Fit Dm Mass to a Non-thermal Profile\n");
	//writeArrayToFile(Rkpc[0], ArrayDividesDouble(dm_mass[0],oneMsunkg), "DMmass");
	dm_density[0] = DensityFromMassSpherical(RLin[0], dm_mass[0]);

	//dm_density[0] = err_one_variable(dm_density[0], dm_mass[0], dm_mass_error[0]);
	//write_array_and_error_to_file(Rkpc[0], dm_density[0], dm_density_error[0], "dMDensity");
	
	mass_jeans[0] = MassFromJeans(RLin[0], dm_density[0], ArraySQ(dm_disp_velocity[0]), dm_beta[0]);
			
	
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(mass_jeans[0],oneMsunkg), "JeansMass");
			
}


		/*array Metropolis_get_bias_jeans_equal_to_he(array radius, array he_mass,  array jeans_mass, array gas_mass, array dm_beta, array dm_disp_vel, 
											array bias_init, int nIter, int nRej, string name_file){
	
	int nBin = radius.numberBins();
	
	array total_mass, dm_density, bias_fit, bias_test, bias_test_bef, bias_test_now;
	total_mass.init(nBin); dm_density.init(nBin); bias_fit.init(nBin); bias_test.init(nBin);bias_test_bef.init(nBin); bias_test_now.init(nBin);
	
	bias_test = CopyArray(bias_init);
	
	double chi_sq[nIter];
	
	//INITIAL CHI_SQ FOR REF
	double chi_sq_test=0.;
	for (int i = 0; i < nBin; i++) {
		double value = pow(he_mass.read(i)-jeans_mass.read(i),2.)/(pow(0.1*he_mass.read(i),2.)+ pow(0.1*jeans_mass.read(i),2.));
		chi_sq_test+= value;
	}
	chi_sq_test = chi_sq_test/nBin;
	printf("chiSQ init = %g \n", chi_sq_test);
	
	chi_sq[0] = chi_sq_test;
	
	// INIT METROPOLIS
	
	double point_a, point_b, chi_ref;
	chi_ref = chi_sq_test;
	
	double tstart, tstop, ttime;
	
	int iter =1; int rej = 0.; int accept = 0;
	int testing; int sum = 0;
	double delta[nBin];
	for (int i = 0; i < nBin; i++) {
		delta[i] = .99;
	}
	ofstream outFile;
	string MyFile = name_file + "_chiSQ.txt";
	outFile.open(MyFile.c_str(), ios::out);
	
	
	while (iter<nIter){
		
		tstart = (double)clock()/CLOCKS_PER_SEC; // starts counting time - we want to know how much time it takes t make the sim
		//printf("METROPOLIS: New configuration %d previous chiSQ = %g when (x0,x1) = (%g,%g) delta = %g\n", iter, chi_sq[iter-1], x0_fit[iter-1], x1_fit[iter-1], delta);
		bool make_new = true;
		while (make_new == true) {
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				bias_test.set(i, bias_test_now.read(i) + cat_rng_uniform_double_min_max(-delta[i],delta[i]));
			}
			total_mass = total_mass_from_bias(he_mass, bias_test);
			dm_density = DensityFromMassSpherical(radius, ArrayMinusArray(total_mass, gas_mass));
			jeans_mass = MassFromJeans(radius, dm_density, ArraySQ(dm_disp_vel), dm_beta);
			
			make_new = false;
			for (int i = 0; i < nBin; i++) {
				if (total_mass.read(i) < gas_mass.read(i)) make_new = true;
				if (dm_density.read(i)!=dm_density.read(i)) make_new = true;
				if (jeans_mass.read(i)!=jeans_mass.read(i)) make_new = true;

			}
			
			for (int j = 0; j < nBin-1; j++) {
				if (total_mass.read(j+1)<=total_mass.read(j)) make_new = true;
			}
			
			for (int j =0 ; j < nBin-1; j++) {
				if (dm_density.read(j)<=dm_density.read(j+1)) make_new = true;
			}
	
			chi_sq_test = 0.;
			for (int j = 0; j < nBin; j++) {
				chi_sq_test+= pow(total_mass.read(j)-jeans_mass.read(j),2.)/(pow(0.1*total_mass.read(j),2.)+ pow(0.1*jeans_mass.read(j),2.));
			}
			chi_sq_test = chi_sq_test/nBin;
			
			if (chi_sq_test != chi_sq_test) make_new = true;
		}
		
	//	printf("chi %g chi bef %g \n", chi_sq_test, chi_sq[iter-1]);
		
		if(chi_sq_test < chi_sq[iter-1]){ // accepts
			chi_sq[iter] = chi_sq_test;
			bias_test_now = CopyArray(bias_test);
			bias_test_bef = CopyArray(bias_test);
			printf("accepting by value: chi_test = %g, chi_before = %g \n", chi_sq_test, chi_sq[iter-1]);
			
			if (chi_sq_test < chi_ref) {
				chi_ref = chi_sq_test;
				bias_fit = CopyArray(bias_test);
			}
			accept++;
		}
		else { //there is a small probability of accepting
			double r = cat_rng_uniform_double_zero_one();
			double w = exp(-.5*(chi_sq_test - chi_sq[iter-1]));
			printf("r = %g   w = %g \n", r, w);
			
			if (r > w) { // i do not accept the new configuration
				chi_sq[iter] = chi_sq[iter-1];
				bias_test_now = CopyArray(bias_test_bef);
				rej++;
				
				if (rej%nRej == 0) delta = 0.9*delta;
			}
			
			if (r<w) {
				chi_sq[iter] = chi_sq_test;
				bias_test_now = CopyArray(bias_test);
				bias_test_bef = CopyArray(bias_test);
				accept++;
			}
			
		}
		
		tstop = (double)clock()/CLOCKS_PER_SEC;
		
		ttime= tstop-tstart; //ttime is how long your code run
		
		printf("It took [%g] time to run one LOOP iter =%d. chi_sq = %g nIter = %d  nAccept = %d\n", ttime, iter, chi_sq_test, nIter, accept);
		
		outFile << iter << " " << chi_sq[iter];
		for (int i = 0; i < nBin; i++) {
			outFile<< " " << bias_test_now.read(i) ;
		}
		 outFile<< endl;
		iter++;
	}
	
	outFile.close();
	
	
	
	
	total_mass.delete_array(); dm_density.delete_array(); 
	
	return bias_fit;
		
	
}
//*/

/*
 Primeiro simples Monte Carlo
 Segundo MonteCarlo com Jump tipo Simmulated Anneling
 Terceiro Genetic Algorithm 
 Por ultimo ver qual deles o melhor e mais rapido
 Há qualquer coisa para contar o tempo que demora - Cpp
 */



