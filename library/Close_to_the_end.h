/*
 *  TheEnd.h
 *  
 *
 *  Created by Catarina Silva Fernandes on 22/01/2011.
 *  Copyright 2011 __MyCompanyName__. All rights reserved.
 *
 */


#include "MyClusterData.h"
/*  
 c++ Close_to_the_end.cpp -o Close_to_the_end.exe -I/local/opt/star/include -L/local/opt/star/lib -l gsl -lgslcblas
 ./Close_to_the_end.exe 
 
 */

//string name_file = argv[1]; 
//double min_plateau = atof(argv[2]);
//double max_plateau = atof(argv[3]);
//double delta = atof(argv[4]);
//int nIter = atoi(argv[5]);
//int nRej = atoi(argv[6]);

array return_deltaM (array he_mass, array bias){
	
	array C; C.init(he_mass.numberBins());
					
	C = ArrayMinusArray(ArrayDividesArray(he_mass, ArrayPlusDouble(bias,1.)), he_mass);
	return C;
}

array save_max_value(array Old, array New){
	int n = Old.numberBins();
	
	array A; A.init(n);
	
	for (int i = 0; i < n; i++) {
		double value;
		if(New.read(i) > Old.read(i)) value = New.read(i);
		if(Old.read(i) > New.read(i)) value = Old.read(i);
		A.set(i, value);
	}
	return A;
	
}
array save_min_value(array Old, array New){
	int n = Old.numberBins();
	
	array A; A.init(n);
	
	for (int i = 0; i < n; i++) {
		
		double value;
		if(New.read(i) < Old.read(i)) value = New.read(i);
		if(Old.read(i) < New.read(i)) value = Old.read(i);
		A.set(i, value);
	}
	return A;
	
}

array give_ref(array Max, array Min){
	int n = Max.numberBins();
	
	array A; A.init(n);
	
	A = ArrayPlusArray(Min, ArrayDividesDouble( ArrayMinusArray( Max,Min),2.));
	
	return A;
}
array give_err(array Max, array Ref){
	int n = Max.numberBins();
	
	array A; A.init(n);
	
	A = ArrayMinusArray(Max, Ref);
	return A;
}
using namespace std;
int main (int argc, char **argv) {
	
	double oneMeterkpc = 1/(3.08578e19); 
	string name_file = argv[1];
	double value_beta = atof(argv[2]);
/*	int nIter = atoi(argv[2]);
	int nRej = atoi(argv[3]);
	double  delta = atof(argv[4]);
	double min_plateau = atof(argv[5]);
	double max_plateau = atof(argv[6]);
//*/	
	//%%%% TAKING A2052 DATA
	ofstream write_file;
	ifstream read_file;
	
	string write_to_file, read_from_file;
	
	srand(time(NULL)); /*genera¥tes seed*/
	
	printf("%g \n", atom*1.e6);
	printf("%g \n", oneMeterkpc);
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//cin >> nIter;
	//string name_file = "t2Wed19_a1689";
	
	write_to_file = name_file;
	take_data_SI_to_files_A2052(write_to_file);
	//	take_cut_data_SI_to_files_A2052(write_to_file);
	//	take_data_SI_to_files_A1689(write_to_file);
	//make fit & print for gas_mass, mass_HE, gas_temperature
	
	
	read_from_file = name_file + "_obs_radius.dat";
	read_file.open(read_from_file.c_str(),ios::in);
	
	double r, r_err; 
	int nBin = 0;
	while (read_file >> r >> r_err) {
		nBin++;
	}
	read_file.close();
	
	array radius_obs, radius_err; radius_obs.init(nBin); radius_err.init(nBin);
	array rkpc; rkpc.init(nBin);
	double r_200 = 1.8e3;
	read_file.open(read_from_file.c_str(),ios::in);
	
	int j = 0;
	while (read_file >> r >> r_err) {
		radius_obs.set(j,r);
		radius_err.set(j,r_err);
		rkpc.set(j, radius_obs.read(j)/( 1.e3*onePcMtr*r_200));
		j++;
	}
	read_file.close();
	
	
	//writeArrayToFile(radius_obs, funct_elephant_inside_snake(radius_obs, radius_obs.read(1), radius_obs.read(2),radius_obs.read(4), radius_obs.read(5), -1., 0.), "_test_snake");
	
	
	//%%%%%%%%%%%%%%%%%%%%%% READ TEMPERATURE DATA %%%%%%%%%%%%%%%%%%%%//
	array gas_temp_obs; gas_temp_obs.init(nBin);
	array gas_temp_err; gas_temp_err.init(nBin);
	
	double par, par_err;
	read_from_file = name_file + "_obs_gas_temp.dat";
	read_file.open(read_from_file.c_str(),ios::in);
	j = 0;
	while (read_file >> r >> par >> par_err) {
		gas_temp_obs.set(j,par);
		gas_temp_err.set(j,par_err);
		j++;
	}
	read_file.close();
	
	//%%%%%%%%%%%%%%%%%%%%%% READ GAS DENSITY DATA %%%%%%%%%%%%%%%%%%%%//
	array gas_dens_obs; gas_dens_obs.init(nBin);
	array gas_dens_obs_err; gas_dens_obs_err.init(nBin);
	
	read_from_file = name_file + "_obs_gas_dens.dat";
	read_file.open(read_from_file.c_str(),ios::in);
	j = 0;
	while (read_file >> r >> par >> par_err) {
		gas_dens_obs.set(j,par);
		gas_dens_obs_err.set(j,par_err);
		j++;
	}
	read_file.close();
	
	
	nBin = 8;
	
	array gas_temp_gauss, gas_dens_gauss, gas_mass_gauss, he_mass_gauss, jeans_mass_gauss;
	gas_temp_gauss.init(nBin); gas_dens_gauss.init(nBin) ; gas_mass_gauss.init(nBin) ; he_mass_gauss.init(nBin)  ; jeans_mass_gauss.init(nBin) ;
	
	array gas_mass_max, gas_mass_min, gas_mass_ref, gas_mass_err, gas_dens_max, gas_dens_min, gas_dens_ref, gas_dens_err;
	gas_mass_max.init(nBin);gas_mass_min.init(nBin);gas_mass_ref.init(nBin);gas_mass_err.init(nBin);
	gas_dens_max.init(nBin); gas_dens_min.init(nBin);gas_dens_ref.init(nBin);gas_dens_err.init(nBin);
	
	
	
	
	array dm_mass_gauss, dm_mass_opt, dm_dens_gauss,dm_dens_opt, dm_beta, dm_disp_vel;
	dm_mass_gauss.init(nBin);dm_dens_gauss.init(nBin); dm_beta.init(nBin) ; dm_disp_vel.init(nBin) ;dm_mass_opt.init(nBin);dm_dens_opt.init(nBin);
	
	array dm_dens_min, dm_dens_max, dm_dens_ref, dm_dens_err;
	dm_dens_min.init(nBin); dm_dens_max.init(nBin); dm_dens_ref.init(nBin); dm_dens_err.init(nBin);
	
	array dm_dens_opt_min, dm_dens_opt_max, dm_dens_opt_ref, dm_dens_opt_err;
	dm_dens_opt_min.init(nBin); dm_dens_opt_max.init(nBin); dm_dens_opt_ref.init(nBin); dm_dens_opt_err.init(nBin);
	
	array dm_mass_ref, dm_mass_max, dm_mass_min, dm_mass_err;
	dm_mass_min.init(nBin); dm_mass_max.init(nBin); dm_mass_ref.init(nBin); dm_mass_err.init(nBin);
	
	array dm_mass_opt_ref, dm_mass_opt_max, dm_mass_opt_min, dm_mass_opt_err;
	dm_mass_opt_ref.init(nBin); dm_mass_opt_max.init(nBin); dm_mass_opt_ref.init(nBin); dm_mass_opt_err.init(nBin);
	
	array jeans_opt, jeans_opt_ref, jeans_opt_min, jeans_opt_max, jeans_opt_err;
	jeans_opt.init(nBin); jeans_opt_ref.init(nBin); jeans_opt_min.init(nBin); jeans_opt_max.init(nBin); jeans_opt_err.init(nBin);
	

	array he_mass_min, he_mass_max, jeans_mass_min, jeans_mass_max,total_mass_min, total_mass_max;
	he_mass_min.init(nBin); he_mass_max.init(nBin); jeans_mass_min.init(nBin); jeans_mass_max.init(nBin); total_mass_min.init(nBin); total_mass_max.init(nBin);
	array total_mass_gauss; total_mass_gauss.init(nBin);
	
	array he_mass_ref, he_mass_err,  total_mass_ref, total_mass_err;
	he_mass_ref.init(nBin); he_mass_err.init(nBin); total_mass_ref.init(nBin); total_mass_err.init(nBin);
	
	array jeans_mass_init, jeans_mass_init_ref, jeans_mass_init_min, jeans_mass_init_max, jeans_mass_init_err;
	jeans_mass_init.init(nBin); jeans_mass_init_ref.init(nBin); jeans_mass_init_min.init(nBin); jeans_mass_init_max.init(nBin);jeans_mass_init_err.init(nBin);	

	
	array bias; bias.init(nBin);
	array bias_ref, bias_err, bias_max, bias_min;
	bias_ref.init(nBin);bias_err.init(nBin);bias_max.init(nBin);bias_min.init(nBin);
	
	//	string MyFile = name_file + "_chi_sq_vs_x0123.txt";
	//	ofstream outFile; outFile.open(MyFile.c_str(), ios::out);
	array bias_err_test, bias_min_test, bias_max_test; bias_err_test.init(nBin); bias_min_test.init(nBin); bias_max_test.init(nBin);
	array total_err_test, total_min_test, total_max_test; total_err_test.init(nBin); total_min_test.init(nBin); total_max_test.init(nBin);
	array dm_err_test, dm_min_test, dm_max_test; dm_err_test.init(nBin); dm_min_test.init(nBin); dm_max_test.init(nBin);
	array jeans_err_test, jeans_min_test, jeans_max_test; jeans_err_test.init(nBin); jeans_min_test.init(nBin); jeans_max_test.init(nBin);
	
	array deltaM, deltaM_ref, deltaM_min, deltaM_max, deltaM_err;
	deltaM.init(nBin); deltaM_ref.init(nBin); deltaM_max.init(nBin); deltaM_err.init(nBin);
	
	int testing;
/*	for(int nGauss = 0; nGauss< 35; nGauss++){
		
		bool make_new = false;
		while (make_new == false) {
		
			for (int j = 0; j < nBin; j ++){
				double value1 = cat_rng_gaussian(1.);
				double value2 = cat_rng_gaussian(1.);
				gas_temp_gauss.set(j, gas_temp_obs.read(j) + value1 * gas_temp_err.read(j));
				gas_dens_gauss.set(j, gas_dens_obs.read(j) + value1 * gas_dens_obs_err.read(j));
			}	
			
			gas_mass_gauss = massFromDensityAssumeSpherical(radius_obs, gas_dens_gauss);
			he_mass_gauss = MassFromHE(radius_obs, gas_dens_gauss, gas_temp_gauss);
			
			dm_mass_gauss = ArrayMinusArray(he_mass_gauss, gas_mass_gauss);
			dm_dens_gauss = DensityFromMassSpherical(radius_obs, dm_mass_gauss);
			
			make_new = true;
			for (int k =0 ; k < nBin; k++) {
				//printf("HE MASS FROM GAUSSIAN %g \n",he_mass_gauss.read(k) );
				if (he_mass_gauss.read(k)<=gas_mass_gauss.read(k)) make_new = false;
			
				if (dm_dens_gauss.read(k)!=dm_dens_gauss.read(k)) make_new = false;
				if (dm_dens_gauss.read(k)==0.) make_new = false;
				
			}
			for (int k =0 ; k < nBin-1; k++) {
			//	if (dm_dens_gauss.read(k)>=dm_dens_gauss.read(k+1)) make_new = false;
			}
			
			int b = -(radius_obs.read(0) + radius_obs.read(1))/(2.*(radius_obs.read(0) - radius_obs.read(1)));
			int m = 1./(radius_obs.read(0) - radius_obs.read(1));
			for (int i = 0; i< nBin; i++) {
			//	dm_beta.set(i,m * radius_obs.read(i) + b);
			}
			//dm_beta = BetaFollowsGammaLinear(Gamma(radius_obs, dm_dens_gauss));
			dm_beta=SetAllToDouble(nBin,value_beta);
			//dm_dens_gauss = DensityFromMassSpherical(radius_obs, ArrayMinusArray(he_mass_gauss, gas_mass_gauss));
			//dm_beta = SetAllToDouble(nBin,0.);	
			
			dm_disp_vel = dispVelFromTemperature(gas_temp_gauss, dm_beta);
			jeans_mass_init = MassFromJeans(radius_obs, dm_dens_gauss, ArraySQ(dm_disp_vel),dm_beta);
			
			for (int k =0 ; k < nBin; k++) {
				if (jeans_mass_init.read(k)!=jeans_mass_init.read(k)) make_new = false;
				if (jeans_mass_init.read(k)<=0.) make_new = false;
			}
			for (int k =0 ; k < nBin-1; k++) {
				if (jeans_mass_init.read(k)>=jeans_mass_init.read(k+1))  make_new = false;
				
			}
		}
		
		
		he_mass_max = save_max_value(he_mass_max, he_mass_gauss); //saves max value (old, new)
		if (nGauss==0) he_mass_min = CopyArray(he_mass_gauss);
		else he_mass_min = save_min_value(he_mass_min, he_mass_gauss);
		
		gas_mass_max =save_max_value(gas_mass_max, gas_mass_gauss); //saves max value (old, new)
		if (nGauss==0) gas_mass_min = CopyArray(gas_mass_gauss);
		else gas_mass_min = save_min_value(gas_mass_min, gas_mass_gauss);
		
		gas_dens_max =save_max_value(gas_dens_max, gas_dens_gauss); //saves max value (old, new)
		if (nGauss==0) gas_dens_min = CopyArray(gas_dens_gauss);
		else gas_dens_min = save_min_value(gas_dens_min,gas_dens_gauss);
		
		jeans_mass_init_max = save_max_value(jeans_mass_init_max, jeans_mass_init); //saves max value (old, new)
		if (nGauss==0) jeans_mass_init_min = CopyArray(jeans_mass_init);
		else jeans_mass_init_min = save_min_value(jeans_mass_init_min, jeans_mass_init);
		
		dm_mass_max = save_max_value(dm_mass_max, dm_mass_gauss); //saves max value (old, new)
		if (nGauss==0) dm_mass_min = CopyArray(dm_mass_gauss);
		else dm_mass_min = save_min_value(dm_mass_min, dm_mass_gauss);
		
		dm_dens_max = save_max_value(dm_dens_max, dm_dens_gauss); //saves max value (old, new)
		if (nGauss==0) dm_dens_min = CopyArray(dm_dens_gauss);
		else dm_dens_min = save_min_value(dm_dens_min, dm_dens_gauss);
		
		
	/*	if (make_new == true) {
			//save previous values
						
			double chi_sq_init = 0.;
			for (int i = 0; i < nBin; i++) {
				chi_sq_init += pow(he_mass_gauss.read(i)-jeans_mass_init.read(i),2.)/(pow(0.1*he_mass_gauss.read(i),2.)+ pow(0.1*jeans_mass_init.read(i),2.));
			}
			
			double * x_left; x_left = new double [2];
			double * x_right; x_right = new double [2];
			
			double * x_left_try; x_left_try = new double [4]; 
			x_left_try[0] = radius_obs.read(0)*0.0001; x_left_try[1] = radius_obs.read(1)*100.;
			x_left_try[2] = min_plateau;x_left_try[3] = max_plateau*0.;
			
			
			double * x_right_try; x_right_try = new double [4]; 
			x_right_try[0] = radius_obs.read(0)*0.01; x_right_try[1] = radius_obs.read(7)*10.;
			x_right_try[2] = min_plateau;x_right_try[3] = max_plateau;
			
			string left = name_file + "_left";
			string right = name_file + "_right";
			
			
			int n_left = 4;
			x_left = A2052_monte_carlo_loop_step_min_max(n_left, radius_obs, he_mass_gauss,jeans_mass_gauss,
														 gas_mass_gauss,  dm_beta,  dm_disp_vel, x_left_try, nIter, nRej, delta,  left);
			
		//	x_left = A2052_monte_carlo_loop_step_max_min_left(n_left, radius_obs, he_mass_gauss,jeans_mass_gauss,
		//												  gas_mass_gauss,  dm_beta,  dm_disp_vel, x_left_try, nIter, nRej, delta,  left);
			int n_right = 4;
			x_right = A2052_monte_carlo_loop_step_max_min(n_right, radius_obs, he_mass_gauss,jeans_mass_gauss,
														  gas_mass_gauss,  dm_beta,  dm_disp_vel, x_right_try, nIter, nRej, delta,  right);
			
			bias = new_funct_biased_elephant_inside_snake (radius_obs, x_left[0], x_left[1], n_left, 
														   x_right[0], x_right[1], n_right, min_plateau, max_plateau);
			
		//	bias = funct_max_to_min(radius_obs, x_right[0], x_right[1], min_plateau, max_plateau);
			
			//bias = metropolis_given_all(radius_obs, he_mass_gauss,  jeans_mass_init,  gas_mass_gauss, dm_beta, dm_disp_vel, 
			//							bias, nIter,  nRej, name_file);
			//bias = Metropolis_get_bias_jeans_equal_to_he(radius_obs, he_mass_gauss,  jeans_mass_init,  gas_mass_gauss, dm_beta, dm_disp_vel, 
			//											 bias, nIter,  nRej, name_file);
			//printf("TEST \n");
			//cin >> testing;
			
			total_mass_gauss = total_mass_from_bias(he_mass_gauss, bias);
			dm_mass_opt = ArrayMinusArray (total_mass_gauss, gas_mass_gauss);
			dm_dens_opt = DensityFromMassSpherical(radius_obs, dm_mass_opt);	
			
			jeans_opt = MassFromJeans(radius_obs, dm_dens_opt, ArraySQ(dm_disp_vel), dm_beta);
			
			double chi_sq_final = 0.;
			for (int i = 0; i < nBin; i++) {
				chi_sq_final += pow(total_mass_gauss.read(i)-jeans_opt.read(i),2.)/(pow(0.1*total_mass_gauss.read(i),2.)+ pow(0.1*jeans_opt.read(i),2.));
				if (jeans_opt.read(i)!= jeans_opt.read(i)) {
					make_new = false;
				}
			}
			if (chi_sq_final != chi_sq_final) {
				make_new = false;
			}
	
		for (int k =0 ; k < nBin-1; k++) {
			if (jeans_opt.read(k)>=jeans_opt.read(k+1))  make_new = false;
			
		}
			
			deltaM = return_deltaM(he_mass_gauss, bias);
			deltaM_max = save_max_value(deltaM_max, deltaM); //saves max value (old, new)
			if (nGauss==0) deltaM_min = CopyArray(deltaM);
			else deltaM_min = save_min_value(deltaM_min, deltaM);
			
			
			total_mass_max = save_max_value(total_mass_max, total_mass_gauss); //saves max value (old, new)
			if (nGauss==0) total_mass_min = CopyArray(total_mass_gauss);
			else total_mass_min = save_min_value(total_mass_min, total_mass_gauss);
			
			
			jeans_opt_max = save_max_value(jeans_opt_max, jeans_opt); //saves max value (old, new)
			if (nGauss==0) jeans_opt_min = CopyArray(jeans_opt);
			else jeans_opt_min = save_min_value(jeans_opt_min, jeans_opt);
			
			dm_mass_opt_max= save_max_value(dm_mass_opt_max, dm_mass_opt); //saves max value (old, new)
			if (nGauss==0) dm_mass_opt_min= CopyArray(dm_mass_opt);
			else dm_mass_opt_min = save_min_value(dm_mass_opt_min, dm_mass_opt);
			
			dm_dens_opt_max = save_max_value(dm_dens_opt_max, dm_dens_opt); //saves max value (old, new)
			if (nGauss==0) dm_dens_opt_min = CopyArray(dm_dens_opt);
			else dm_dens_opt_min = save_min_value(dm_dens_opt_min, dm_dens_opt);
			
			bias_max = save_max_value(bias_max, bias); //saves max value (old, new)
			if (nGauss==0) bias_min = CopyArray(bias);
			else bias_min = save_min_value(bias_min, bias);
			}
	 //*/
	 //*/}
	
	printf("HEJ\n");
	
	he_mass_ref = give_ref (he_mass_max, he_mass_min);
	he_mass_err = give_err (he_mass_max, he_mass_ref);
	
	gas_mass_ref = give_ref (gas_mass_max, gas_mass_min);
	gas_mass_err = give_err (gas_mass_max, gas_mass_ref);

	
	jeans_mass_init_ref = give_ref (jeans_mass_init_max, jeans_mass_init_min);
	jeans_mass_init_err = give_err (jeans_mass_init_max, jeans_mass_init_ref);
	
	
	dm_mass_ref = give_ref (dm_mass_max, dm_mass_min);
	dm_mass_err = give_err (dm_mass_max, dm_mass_ref);
	
//	dm_mass_opt_ref = give_ref (dm_mass_opt_max, dm_mass_opt_min);
//	dm_dens_opt_err = give_err (dm_mass_opt_max, dm_mass_opt_ref);
	
	//fit to DM density
//*/	
	dm_dens_ref = give_ref( dm_dens_max, dm_dens_min);
	dm_dens_err = give_err (dm_dens_max, dm_dens_ref);

	dm_dens_opt_ref = give_ref( dm_dens_opt_max, dm_dens_opt_min);
	dm_dens_opt_err = give_err (dm_dens_opt_max, dm_dens_opt_ref);
	
	
	printf("hejj \n");
/*	double rho0 = 1.2e-22  ; double r_s = 2.0e21; double alpha = .5; double beta = 3.0;
	double dm_dens_guess[4] = { 10.*rho0, r_s, alpha, beta};
	double * fit_parameters; fit_parameters = new double [8];
	fit_parameters = fit_density_PL_w_error_return_parameters_set_alpha_beta(radius_obs, dm_dens_ref, dm_dens_err, dm_dens_guess);
	
	printf("DM init parameters = (%g, %g) power (%g, %g)\n", fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	printf("DM error = (%g, %g) power (%g, %g)\n", fit_parameters[4], fit_parameters[5], fit_parameters[6], fit_parameters[7], fit_parameters[8]);
	
	fit_parameters = fit_density_PL_w_error_return_parameters_set_alpha_beta(radius_obs, dm_dens_opt_ref, dm_dens_opt_err, dm_dens_guess);
	printf("DM final parameters = (%g, %g) power (%g, %g)\n", fit_parameters[0], fit_parameters[1], fit_parameters[2], fit_parameters[3]);
	printf("DM error = (%g, %g) power (%g, %g)\n", fit_parameters[4], fit_parameters[5], fit_parameters[6], fit_parameters[7], fit_parameters[8]);
	
//*/	
	printf("hejj \n");
	jeans_mass_init_ref = give_ref(jeans_mass_init_max, jeans_mass_init_min);
	jeans_mass_init_err = give_err(jeans_mass_init_max, jeans_mass_init_ref);
	ofstream outFile_he, outFile_gas_mass, outFile_jeans_init, outFile_jeans_opt, outFile_dm_mass_init, outFile_total, outFile_deltaM;
	
	
	string HE_name = name_file+"_he.txt";
	outFile_he.open(HE_name.c_str(), ios::out);
	
	string gas_mass_init_name = name_file+"_gas_mass.txt";
	outFile_gas_mass.open(gas_mass_init_name.c_str(), ios::out);
	
	
	string dm_mass_init_name = name_file+"_dm_mass_init.txt";
	outFile_dm_mass_init.open(dm_mass_init_name.c_str(), ios::out);
	
	string Jeans_init_name = name_file+"_jeans_init.txt";
	outFile_jeans_init.open(Jeans_init_name.c_str(), ios::out);
	
	
	// save mass
	for (int j = 0; j < nBin; j++) {
		outFile_he << radius_obs.read(j)*oneMeterkpc << " " << he_mass_ref.read(j)/oneMsunkg << " " << he_mass_err.read(j)/oneMsunkg<< endl;
		outFile_jeans_init << radius_obs.read(j)*oneMeterkpc << " " << jeans_mass_init_ref.read(j)/oneMsunkg << " " << jeans_mass_init_err.read(j)/oneMsunkg<< endl;
		outFile_dm_mass_init << radius_obs.read(j)*oneMeterkpc << " " << dm_mass_ref.read(j)/oneMsunkg << " " << dm_mass_err.read(j)/oneMsunkg<< endl;
		outFile_gas_mass << radius_obs.read(j)*oneMeterkpc << " " << gas_mass_ref.read(j)/oneMsunkg << " " << gas_mass_err.read(j)/oneMsunkg<< endl;
	}
	outFile_he.close();
	outFile_jeans_init.close();
	outFile_dm_mass_init.close();
	outFile_gas_mass.close();
	
	
	double chi_i=0.;
	for (int i = 0; i < nBin; i++) {
		//chi_f += pow(total_mass_ref.read(i)-jeans_opt_ref.read(i),2.)/(pow(total_mass_err.read(i),2.)+ pow(jeans_opt_err.read(i),2.));
		chi_i += pow(he_mass_ref.read(i)-jeans_mass_init_ref.read(i),2.)/(pow(he_mass_err.read(i),2.)+ pow(jeans_mass_init_err.read(i),2.));
		
	}
	
	printf("chi = %g \n", chi_i);
	
	
	
	
	ifstream inpFile, inpFile2;
	inpFile.open("new_betaZero_Man24_jeans_mass_err.txt", ios::in);
	
	array raio, bias_0, bias_0_err; raio.init(nBin); bias_0.init(nBin); bias_0_err.init(nBin);
	array total_0, total_0_err; total_0.init(nBin); total_0_err.init(nBin);
	
	double rad, val, val_err;
	int i =0;
	
	while (inpFile >> rad>> val>> val_err) {
		bias_0.set(i, val);
		bias_0_err.set(i, val_err);
		i++;
	}
	inpFile.close();
	inpFile2.open("new_betaZero_Man24_total_mass_err.txt", ios::in);
	while (inpFile2 >> rad>> val>> val_err) {
		total_0.set(i, val);
		total_0_err.set(i, val_err);
		i++;
	}
	inpFile2.close();
	
	inpFile.open("new_beta13_Man24_jeans_mass_err.txt", ios::in);
	array bias_13, bias_13_err; bias_13.init(nBin); bias_13_err.init(nBin);
	array total_13, total_13_err; total_13.init(nBin); total_13_err.init(nBin);
	
	i =0;
	while (inpFile >> rad>> val>> val_err) {
		
		raio.set(i, rad);
		bias_13.set(i, val);
		bias_13_err.set(i, val_err);
		i++;
	}
	inpFile.close();
	inpFile2.open("new_beta13_Man24_total_mass_err.txt", ios::in);
	
	while (inpFile2 >> rad>> val>> val_err) {
		
		raio.set(i, rad);
		total_13.set(i, val);
		total_13_err.set(i, val_err);
		i++;
	}
	inpFile2.close();
	
	array bias_tru, bias_tru_err; bias_tru.init(nBin); bias_tru_err.init(nBin);
	array total_tru, total_tru_err; total_tru.init(nBin); total_tru_err.init(nBin);
	
	//bias_tru_err = ArraySQRT(ArrayPlusArray(ArraySQ(save_max_value(bias_0_err, bias_13_err)),
	//										ArraySQ(ArrayAbs(DoubleTimesArray(2.,ArrayMinusArray(bias_0,bias_13))) )));
	bias_tru_err = ArraySQRT(ArrayPlusArray(ArraySQ(save_max_value(bias_0_err, bias_13_err)),
											ArraySQ(ArrayAbs(ArrayMinusArray(bias_13,bias_0)))));
	
	bias_tru = ArrayPlusArray(bias_0, ArrayAbs(ArrayDividesDouble(ArrayMinusArray(bias_13, bias_0),2.)));
	
	bias_tru_err = ArrayMinusArray(bias_tru_err, bias_tru);
	bias_tru_err = CopyArray(bias_0_err);
	
	total_tru_err = ArraySQRT(ArrayPlusArray(ArraySQ(ArraySQRT(ArrayPlusArray(ArraySQ(total_0_err), ArraySQ(total_13_err)))),
											ArraySQ(ArrayAbs(DoubleTimesArray(1.,ArrayMinusArray(total_0,total_13))) )));
	
	total_tru = ArrayPlusArray(bias_0, ArrayDividesDouble(ArrayAbs(ArrayMinusArray(total_13, total_0)),2.));
	
	ofstream outFile_bias_get_it, outFile_total_get_it;
	string bias_get_it = name_file+"_jeans_revisited.txt";
	string total_get_it = name_file+"_total_revisited.txt";
	
	outFile_bias_get_it.open(bias_get_it.c_str(),ios::out);
	outFile_total_get_it.open(total_get_it.c_str(),ios::out);
	for (int j = 0; j < nBin; j++) {
		outFile_bias_get_it << raio.read(j) << " " << bias_tru.read(j) << " " << bias_tru_err.read(j)<< endl;
		outFile_total_get_it << raio.read(j) << " " << total_tru.read(j) << " " << total_tru_err.read(j)<< endl;
	}	

	outFile_bias_get_it.close();
/*	jeans_opt_ref = give_ref (jeans_opt_max, jeans_opt_min);
	jeans_opt_err = give_err (jeans_opt_max, jeans_opt_ref);
	
	total_mass_ref = give_ref (total_mass_max, total_mass_min);
	total_mass_err = give_err (total_mass_max, total_mass_ref);
	
	bias_ref = give_ref (bias_max, bias_min);
	bias_err = give_err (bias_max, bias_ref);
	
	deltaM_ref =  give_ref (deltaM_max, deltaM_min);
	deltaM_err =  give_err (deltaM_max, deltaM_ref);
	
	ofstream outFile_he, outFile_gas_mass, outFile_jeans_init, outFile_jeans_opt, outFile_dm_mass_init, outFile_total, outFile_deltaM;
	
	
	
	string Jeans_opt_name = name_file+"_jeans_opt.txt";
	outFile_jeans_opt.open(Jeans_opt_name.c_str(), ios::out);

	string total_name = name_file+"_total_opt.txt";
	outFile_total.open( total_name.c_str(), ios::out);
	string deltaM_name = name_file+"_deltaM.txt";
	outFile_deltaM.open(deltaM_name.c_str(), ios::out);
	
	
	
	
	// save mass
	for (int j = 0; j < nBin; j++) {
		outFile_he << radius_obs.read(j)*oneMeterkpc << " " << he_mass_ref.read(j)/oneMsunkg << " " << he_mass_err.read(j)/oneMsunkg<< endl;
		outFile_jeans_init << radius_obs.read(j)*oneMeterkpc << " " << jeans_mass_init_ref.read(j)/oneMsunkg << " " << jeans_mass_init_err.read(j)/oneMsunkg<< endl;
		outFile_jeans_opt << radius_obs.read(j)*oneMeterkpc << " " << jeans_opt_ref.read(j)/oneMsunkg << " " << jeans_opt_err.read(j)/oneMsunkg<< endl;
		outFile_dm_mass_init << radius_obs.read(j)*oneMeterkpc << " " << dm_mass_ref.read(j)/oneMsunkg << " " << dm_mass_err.read(j)/oneMsunkg<< endl;
		outFile_gas_mass << radius_obs.read(j)*oneMeterkpc << " " << gas_mass_ref.read(j)/oneMsunkg << " " << gas_mass_err.read(j)/oneMsunkg<< endl;
		outFile_total << radius_obs.read(j)*oneMeterkpc << " " << total_mass_ref.read(j)/oneMsunkg << " " << total_mass_err.read(j)/oneMsunkg<< endl;
		outFile_deltaM << radius_obs.read(j)*oneMeterkpc << " " << deltaM_ref.read(j)/oneMsunkg << " " << deltaM_err.read(j)/oneMsunkg<< endl;
		
	}
	outFile_he.close();
	outFile_jeans_init.close();
	outFile_jeans_opt.close();
	outFile_dm_mass_init.close();
	outFile_gas_mass.close();
	outFile_total.close();
	outFile_deltaM.close();
	
	
	ofstream outFile_dm_dens_init, outFile_gas_dens, outFile_dm_dens_opt, outFile_bias; 
	
	string dm_dens_init_name = name_file+"_dm_dens_init.txt";
	outFile_dm_dens_init.open(dm_dens_init_name.c_str(), ios::out);
	
	string gas_dens_name = name_file+"_gas_dens.txt";
	outFile_gas_dens.open(gas_dens_name.c_str(), ios::out);
	
	string dm_dens_opt_name = name_file+"_dm_dens_opt.txt";
	outFile_dm_dens_opt.open(dm_dens_opt_name.c_str(), ios::out);
	
	string bias_name = name_file+"_bias.txt";
	outFile_bias.open(bias_name.c_str(), ios::out);
	

	// save mass
	for (int j = 0; j < nBin; j++) {
		outFile_dm_dens_init << radius_obs.read(j)*oneMeterkpc << " " << dm_dens_ref.read(j)/(atom*1.e6) << " " << dm_dens_err.read(j)/(atom*1.e6)<< endl;
		outFile_gas_dens << radius_obs.read(j)*oneMeterkpc << " " << gas_dens_ref.read(j)/(atom*1.e6) << " " << gas_dens_err.read(j)/(atom*1.e6)<< endl;
		outFile_dm_dens_opt << radius_obs.read(j)*oneMeterkpc << " " << dm_dens_opt_ref.read(j)/(atom*1.e6)<< " " << dm_dens_opt_err.read(j)/(atom*1.e6)<< endl;
		outFile_bias << radius_obs.read(j)*oneMeterkpc << " " << bias_ref.read(j)*100. << " " << bias_err.read(j)*100. << endl;
	}
	outFile_dm_dens_init.close();
	outFile_gas_dens.close();
	outFile_dm_dens_opt.close();
	outFile_bias.close();
	
	
	double chi_f=0., chi_i=0.;
	for (int i = 0; i < nBin; i++) {
		chi_f += pow(total_mass_ref.read(i)-jeans_opt_ref.read(i),2.)/(pow(total_mass_err.read(i),2.)+ pow(jeans_opt_err.read(i),2.));
		chi_i += pow(he_mass_ref.read(i)-jeans_mass_init_ref.read(i),2.)/(pow(he_mass_err.read(i),2.)+ pow(jeans_mass_init_err.read(i),2.));
		
	}
	
	
	printf("init %g final  %g \n", chi_i, chi_f);
	//*/
	
	
	inpFile.open("new_DMREADMEBetaZero.txt", ios::in);
	i=0;
	while (inpFile >> rad>> val>> val_err) {
		i++;
	}
	inpFile.close();
	
	array R_gas, Gas, gamma;
	R_gas.init(i); Gas.init(i); gamma.init(i);
	i =0;
	
	inpFile.open("new_DMREADMEBetaZero.txt", ios::in);
	while (inpFile >> rad>> val>> val_err) {
		R_gas.set(i, rad);
		Gas.set(i, val);
		i++;
	}
	inpFile.close();
	
	gamma = Gamma(R_gas, Gas);
	
	writeArrayToFile(R_gas, gamma, "Gamma");
	
	
	return 0;

}