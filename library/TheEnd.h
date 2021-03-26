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
 c++ TheEnd.cpp -o TheEnd.exe -I/local/opt/star/include -L/local/opt/star/lib -l gsl -lgslcblas
 ./TheEnd.exe name_file
 */

//string name_file = argv[1]; 
//double min_plateau = atof(argv[2]);
//double max_plateau = atof(argv[3]);
//double delta = atof(argv[4]);
//int nIter = atoi(argv[5]);
//int nRej = atoi(argv[6]);


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
	//cout<< " hhhh " << atom << endl;
	// Choose number of Interations, Rejections, init Delta, x{2}, plateau{}
	//int nIter = atoi(argv[1]);
	//int nRej = atoi(argv[2]);
	//double  delta = atof(argv[3]);
	//double x_i[2] = {atof(argv[4]), atof(argv[5])};
	string name_file = argv[1];
	
	//%%%% TAKING A2052 DATA
	ofstream write_file;
	ifstream read_file;
	
	string write_to_file, read_from_file;
	
	srand(time(NULL)); /*genera¥tes seed*/
	
	
	
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
	array gas_dens_err; gas_dens_err.init(nBin);
	
	read_from_file = name_file + "_obs_gas_dens.dat";
	read_file.open(read_from_file.c_str(),ios::in);
	j = 0;
	while (read_file >> r >> par >> par_err) {
		gas_dens_obs.set(j,par);
		gas_dens_err.set(j,par_err);
		j++;
	}
	read_file.close();
	
	
	nBin = 8;
	
	array gas_temp_gauss, gas_dens_gauss, gas_mass_gauss, he_mass_gauss, jeans_mass_gauss;
	gas_temp_gauss.init(nBin); gas_dens_gauss.init(nBin) ; gas_mass_gauss.init(nBin) ; he_mass_gauss.init(nBin)  ; jeans_mass_gauss.init(nBin) ;
	
	array dm_dens_gauss, dm_beta, dm_disp_vel;
	dm_dens_gauss.init(nBin); dm_beta.init(nBin) ; dm_disp_vel.init(nBin) ;
	
	dm_beta = SetAllToDouble(nBin,0.);	
	
	ofstream outFile0, outFile1, outFile2, outFile3, outFile4;
	string name_he = name_file+"he_gauss.txt";
	string name_jeans = name_file+"jeans_gauss.txt";
	string name_gas_mass = name_file+"gas_mass.txt";
	string name_dm_dens = name_file+"dm_dens_gauss.txt";
	
//	outFile0.open(name_he.c_str(), ios::out);
//	outFile1.open(name_jeans.c_str(), ios::out);
//	outFile2.open(name_gas_mass.c_str(), ios::out);
	//outFile3.open(name_dm_dens.c_str(), ios::out);
	
	array jeans_mass_min_out, jeans_mass_max_out, jeans_mass_ref_out;
	array  total_mass_max_out,  total_mass_min_out, total_mass_ref_out;
	
	total_mass_ref_out.init(nBin); total_mass_max_out.init(nBin); total_mass_min_out.init(nBin); 
	jeans_mass_ref_out.init(nBin); jeans_mass_max_out.init(nBin); jeans_mass_min_out.init(nBin); 
	
	array jeans_mass_err_out, total_mass_err_out;
	jeans_mass_err_out.init(nBin); total_mass_err_out.init(nBin);
	
	array he_mass_min, he_mass_max, jeans_mass_min, jeans_mass_max,total_mass_min, total_mass_max;
	he_mass_min.init(nBin); he_mass_max.init(nBin); jeans_mass_min.init(nBin); jeans_mass_max.init(nBin); total_mass_min.init(nBin); total_mass_max.init(nBin);
	array total_mass_gauss; total_mass_gauss.init(nBin);
	
	array dm_dens_min, dm_dens_max, dm_dens_ref, dm_dens_err;
	dm_dens_min.init(nBin); dm_dens_max.init(nBin); dm_dens_ref.init(nBin); dm_dens_err.init(nBin);
	
	array he_mass_ref, he_mass_err, jeans_mass_ref, jeans_mass_err, total_mass_ref, total_mass_err;
	he_mass_ref.init(nBin); he_mass_err.init(nBin); jeans_mass_ref.init(nBin); jeans_mass_err.init(nBin);total_mass_ref.init(nBin); total_mass_err.init(nBin);
	
	array bias; bias.init(nBin);
	array bias_ref, bias_err, bias_max, bias_min;
	bias_ref.init(nBin);bias_err.init(nBin);bias_max.init(nBin);bias_min.init(nBin);
	
//	string MyFile = name_file + "_chi_sq_vs_x0123.txt";
//	ofstream outFile; outFile.open(MyFile.c_str(), ios::out);
	array bias_err_test, bias_min_test, bias_max_test; bias_err_test.init(nBin); bias_min_test.init(nBin); bias_max_test.init(nBin);
	array total_err_test, total_min_test, total_max_test; total_err_test.init(nBin); total_min_test.init(nBin); total_max_test.init(nBin);
	array dm_err_test, dm_min_test, dm_max_test; dm_err_test.init(nBin); dm_min_test.init(nBin); dm_max_test.init(nBin);
	array jeans_err_test, jeans_min_test, jeans_max_test; jeans_err_test.init(nBin); jeans_min_test.init(nBin); jeans_max_test.init(nBin);
	array gas_mass_err; gas_mass_err.init(nBin);
	int testing;
	for(int nGauss = 0; nGauss< 35; nGauss++){
		
	
		
		bool make_new = false;
	/*	while (make_new == false) {
			for (int j = 0; j < nBin; j ++){
				double value1 = cat_rng_gaussian(1.);
				double value2 = cat_rng_gaussian(1.);
				gas_temp_gauss.set(j, gas_temp_obs.read(j) + value1 * gas_temp_err.read(j));
				gas_dens_gauss.set(j, gas_dens_obs.read(j) + value1 * gas_dens_err.read(j));
			}	
				
			gas_mass_gauss = massFromDensityAssumeSpherical(radius_obs, gas_dens_gauss);
			he_mass_gauss = MassFromHE(radius_obs, gas_dens_gauss, gas_temp_gauss);
			
			make_new = true;
			for (int k =0 ; k < nBin; k++) {
				//printf("HE MASS FROM GAUSSIAN %g \n",he_mass_gauss.read(k) );
				if (he_mass_gauss.read(k)<=gas_mass_gauss.read(k)) {
					make_new = false;
				}
				if ( DensityFromMassSpherical(radius_obs,he_mass_gauss).read(k)<=gas_dens_gauss.read(k)) {
					make_new = false;
				}
				if (DensityFromMassSpherical(radius_obs,he_mass_gauss).read(k)<0.) {
					make_new = false;
				}
				printf("HE dens = %g   Gas dens =%g\n", DensityFromMassSpherical(radius_obs,he_mass_gauss).read(k), gas_dens_gauss.read(k));
			}
			//dm_beta = BetaFollowsGammaLinear(Gamma(radius_obs, dm_dens_gauss));
			//dm_beta=SetAllToDouble(nBin,1./3.);
			//dm_dens_gauss = DensityFromMassSpherical(radius_obs, ArrayMinusArray(he_mass_gauss, gas_mass_gauss));
			
			
			
			dm_dens_gauss = ArrayMinusArray(he_mass_gauss, gas_mass_gauss);
			for (int k =0 ; k < nBin-1; k++) {
				if (dm_dens_gauss.read(k)>=dm_dens_gauss.read(k+1)) {
					make_new = false;
				}
				
			}
			for (int k =0 ; k < nBin; k++) {
				if (dm_dens_gauss.read(k)!=dm_dens_gauss.read(k)) {
					make_new = false;
				}
				if (dm_dens_gauss.read(k)<gas_mass_gauss.read(k)) {
					make_new = false;
				}
			}
		/*	for (int k =0 ; k < nBin; k++) {
				if (dm_dens_gauss.read(k)<gas_dens_gauss.read(k)) {
					make_new = false;
					printf("k=%d\n",k);
				}
			}
		/*	dm_disp_vel = dispVelFromTemperature(gas_temp_gauss, dm_beta);
			
			jeans_mass_gauss = MassFromJeans(radius_obs, dm_dens_gauss, ArraySQ(dm_disp_vel),dm_beta);
			for (int k =0 ; k < nBin; k++) {
				if (jeans_mass_gauss.read(k)!=jeans_mass_gauss.read(k)) {
					make_new = false;
				}
			}
			for (int k =0 ; k < nBin-1; k++) {
				if (jeans_mass_gauss.read(k)>=jeans_mass_gauss.read(k+1)) {
					make_new = false;
				}
				if (jeans_mass_gauss.read(k)<0.) {
					make_new = false;
				}
			}
		//*/	
		/*	if (make_new == true) {
				he_mass_max = save_max_value(he_mass_max, he_mass_gauss); //saves max value (old, new)
				if (nGauss==0) he_mass_min = CopyArray(he_mass_gauss);
				else he_mass_min = save_min_value(he_mass_min, he_mass_gauss);
				
				bias_min = save_min_value(bias_min, bias);
				if (nGauss==0)bias_max = CopyArray(bias);
				else bias_max = save_max_value(bias_max, bias);
				
				jeans_mass_max = save_max_value(jeans_mass_max, jeans_mass_gauss);
				if (nGauss==0) jeans_mass_min = CopyArray(jeans_mass_gauss);
				else jeans_mass_min = save_min_value(jeans_mass_min, jeans_mass_gauss);
				
			/*	total_mass_max = save_max_value(total_mass_max, total_mass_gauss);
				if (nGauss==0) total_mass_min = CopyArray(total_mass_gauss);
				else total_mass_min = save_min_value(total_mass_min, total_mass_gauss);
			//*/	
			/*	dm_dens_max = save_max_value(dm_dens_max, dm_dens_gauss);
				if (nGauss==0) dm_dens_min = CopyArray(dm_dens_gauss);
				else dm_dens_min = save_min_value(dm_dens_min, dm_dens_gauss);
				
				total_mass_max = save_max_value(total_mass_max, gas_mass_gauss);
				if (nGauss==0) total_mass_min = CopyArray(gas_mass_gauss);
				else total_mass_min = save_min_value(total_mass_min, gas_mass_gauss);
				
			}
		//*/	
		/*	if (make_new == true) {
				double chi_sq_init = 0.;
				for (int i = 0; i < nBin; i++) {
					chi_sq_init += pow(he_mass_gauss.read(i)-jeans_mass_gauss.read(i),2.)/(pow(0.1*he_mass_gauss.read(i),2.)+ pow(0.1*jeans_mass_gauss.read(i),2.));
				}
				chi_sq_init = chi_sq_init;
				
				
				double min_plateau = atof(argv[2]); 
				double max_plateau = atof(argv[3]);
				//double x_guess[6] = {1.e20, radius_obs.read(4), radius_obs.read(5), 1.e23,min_plateau, max_plateau};
				//double x_guess[6] = {1.e20, radius_obs.read(4), radius_obs.read(5),  radius_obs.read(7),min_plateau, max_plateau};
				double x_guess[6] = {radius_obs.read(0), radius_obs.read(3), radius_obs.read(5),  radius_obs.read(7),min_plateau, max_plateau};
				double * x_left_try; x_left_try = new double [4]; 
				x_left_try[0] = x_guess[0]; x_left_try[1] = x_guess[1];x_left_try[2] = x_guess[4];x_left_try[3] = x_guess[5];
				double * x_right_try; x_right_try = new double [4];
				x_right_try[0] = x_guess[2]; x_right_try[1] = x_guess[3]; x_right_try[2] = x_guess[4];x_right_try[3] = x_guess[5];
				
				double * x_left; x_left = new double [2];
				double * x_right; x_right = new double [2];
				int nIter = atoi(argv[5]); int nRej = atoi(argv[6]); double delta = atof(argv[4]);
				
				//x_best = Metropolis_get_jeans_equal_to_he( radius_obs,  he_mass_gauss,  he_mass_err,  jeans_mass_gauss,  jeans_mass_err,
				//										  gas_mass_gauss,  dm_beta,  dm_disp_vel, x_guess,  nIter,  nRej, delta, name_file);
				
				int n_left = 5;
				x_left = A2052_monte_carlo_loop_step_min_max(n_left, radius_obs, he_mass_gauss,jeans_mass_gauss,
															 gas_mass_gauss,  dm_beta,  dm_disp_vel, x_left_try, nIter, nRej, delta,  name_file);
				
				 min_plateau = atof(argv[7]);
				 max_plateau = atof(argv[8]);
				
				int n_right = 5.;
				x_right = A2052_monte_carlo_loop_step_max_min(n_right, radius_obs, he_mass_gauss,jeans_mass_gauss,
					  gas_mass_gauss,  dm_beta,  dm_disp_vel, x_right_try, nIter, nRej, delta,  name_file);
				
				//bias = funct_min_to_max(radius_obs, x_left[0], x_left[1], min_plateau, max_plateau);
				//bias = funct_elephant_inside_snake(radius_obs, x_left[0], x_left[1], x_right[0], x_right[1], min_plateau, max_plateau);
				bias = new_funct_biased_elephant_inside_snake (radius_obs, x_left[0], x_left[1], n_left, 
															   x_right[0], x_right[1], n_right, min_plateau, max_plateau);
				
				bias_err_test = DoubleTimesArray(0.1, bias);
				
				bias_min_test = ArrayMinusArray(bias, bias_err_test);
				bias_max_test = ArrayPlusArray(bias, bias_err_test);
			//*/	
				
		/*		total_mass_gauss = total_mass_from_bias(he_mass_gauss, bias);
				
				total_err_test = ArraySQRT(ArraySQ(ArrayTimesArray(Derivative(bias, total_mass_gauss), bias_err_test)));
				total_min_test= ArrayMinusArray(total_mass_gauss, total_err_test);
				total_max_test = ArrayPlusArray(total_mass_gauss, total_err_test);

				
				
				dm_dens_gauss = DensityFromMassSpherical(radius_obs, ArrayMinusArray (total_mass_gauss, gas_mass_gauss));		
				dm_err_test = ArraySQRT(ArraySQ(ArrayTimesArray(Derivative(total_mass_gauss, dm_dens_gauss), total_err_test)));
				dm_min_test= ArrayMinusArray(dm_dens_gauss, dm_err_test);
				dm_max_test = ArrayPlusArray(dm_dens_gauss, dm_err_test);
							
				
				jeans_mass_gauss = MassFromJeans(radius_obs, dm_dens_gauss, ArraySQ(dm_disp_vel), dm_beta);
				jeans_err_test = ArraySQRT(ArraySQ(ArrayTimesArray(Derivative(dm_dens_gauss, jeans_mass_gauss), dm_err_test)));
				jeans_min_test= ArrayMinusArray(jeans_mass_gauss, jeans_err_test);
				jeans_max_test = ArrayPlusArray(jeans_mass_gauss, jeans_err_test);
				
	
				
				for (int k = 0; k < nBin; k++) {
					if (total_mass_gauss.read(k)!=total_mass_gauss.read(k)) make_new = false;
					if (jeans_mass_gauss.read(k)!=jeans_mass_gauss.read(k)) make_new = false;
					if (dm_dens_gauss.read(k)!=dm_dens_gauss.read(k)) make_new = false;
					
					if (k < nBin-1) {
						if (total_mass_gauss.read(k) > total_mass_gauss.read(k+1)) make_new = false;
						if (jeans_mass_gauss.read(k) > jeans_mass_gauss.read(k+1)) make_new = false;
						if (dm_dens_gauss.read(k) < dm_dens_gauss.read(k+1)) make_new = false;
					}
				}
				
				
				double chi_sq_final = 0.;
				for (int i = 0; i < nBin; i++) {
					chi_sq_final += pow(total_mass_gauss.read(i)-jeans_mass_gauss.read(i),2.)/(pow(total_err_test.read(i),2.)+ pow(jeans_err_test.read(i),2.));
				}
				chi_sq_final = chi_sq_final;
				
				if (chi_sq_init <= chi_sq_final || chi_sq_final != chi_sq_final) make_new = false;
				printf("chi_init = %g , chi_final = %g \n", chi_sq_init, chi_sq_final);
				
				if (make_new == true) {
					he_mass_max = save_max_value(he_mass_max, he_mass_gauss); //saves max value (old, new)
					if (nGauss==0) he_mass_min = CopyArray(he_mass_gauss);
					else he_mass_min = save_min_value(he_mass_min, he_mass_gauss);

				/*	bias_min = save_min_value(bias_min, bias);
					if (nGauss==0)bias_max = CopyArray(bias);
					else bias_max = save_max_value(bias_max, bias);
				//*/
					
				/*	jeans_mass_max = save_max_value(jeans_mass_max, jeans_mass_gauss);
					if (nGauss==0) jeans_mass_min = CopyArray(jeans_mass_gauss);
					else jeans_mass_min = save_min_value(jeans_mass_min, jeans_mass_gauss);
				//*/	
			//		total_mass_max = save_max_value(total_mass_max, gas_dens_gauss);
			//		if (nGauss==0) total_mass_min = CopyArray(gas_dens_gauss);
			//		else total_mass_min = save_min_value(total_mass_min, gas_dens_gauss);
				//*/	
			//		dm_dens_max = save_max_value(dm_dens_max, dm_dens_gauss);
			//		if (nGauss==0) dm_dens_min = CopyArray(dm_dens_gauss);
			//		else dm_dens_min = save_min_value(dm_dens_min, dm_dens_gauss);
					
				//*/
					
				/*	jeans_mass_max = save_max_value(jeans_mass_max, jeans_max_test);
					if (nGauss==0) jeans_mass_min = CopyArray(jeans_mass_gauss);
					else jeans_mass_min = save_min_value(jeans_mass_min, jeans_min_test);
					
					total_mass_max = save_max_value(total_mass_max, total_max_test);
					if (nGauss==0) total_mass_min = CopyArray(total_mass_gauss);
					else total_mass_min = save_min_value(total_mass_min, total_min_test);
					
					dm_dens_max = save_max_value(dm_dens_max, dm_max_test);
					if (nGauss==0) dm_dens_min = CopyArray(dm_dens_gauss);
					else dm_dens_min = save_min_value(dm_dens_min, dm_min_test);
					
					
					
					bias_min = save_min_value(bias_min, bias_min_test);
					if (nGauss==0)bias_max = CopyArray(bias);
					else bias_max = save_max_value(bias_max, bias_max_test);
					//*/
					
					
					//cin >> testing;
					
				}
			
		//*/	
	//*/	}
		//gas dens
	//	jeans_mass_max = save_max_value(jeans_mass_max, gas_dens_gauss);
	//	if (nGauss==0) jeans_mass_min = CopyArray(gas_dens_gauss);
	//	else jeans_mass_min = save_min_value(jeans_mass_min, gas_dens_gauss);
		
		//dm dens
	//	total_mass_max = save_max_value(total_mass_max, dm_dens_gauss);
	//	if (nGauss==0) total_mass_min = CopyArray(dm_dens_gauss);
	//	else total_mass_min = save_min_value(total_mass_min, dm_dens_gauss);
		
		
		
		
		//bias =  funct_elephant_inside_snake(radius_obs, x_left[0],x_left[1],x_right[2],x_right[3], min_plateau, max_plateau);
		 
		
		//writeArrayToFile(radius_obs, bias, "_bias");
	//	total_mass_gauss = total_mass_from_bias(he_mass_gauss, bias);
		
	//	
		

		
		//for (int j = 0; j < nBin; j++) {
			//outFile0 << radius_obs.read(j) << " " << he_mass_gauss.read(j) << endl;
			//outFile1 << radius_obs.read(j) << " " << jeans_mass_gauss.read(j) << endl;
			//outFile2 << radius_obs.read(j) << " " << total_mass_gauss.read(j) << endl;
		//}
		//*/
	//}
	//outFile.close();
	//outFile0.close();
	//outFile1.close();
	//outFile2.close();
	
	
	he_mass_ref = give_ref (he_mass_max, he_mass_min);
	he_mass_err = give_err (he_mass_max, he_mass_ref);
	
	jeans_mass_ref = give_ref (jeans_mass_max, jeans_mass_min);
	jeans_mass_err = give_err (jeans_mass_max, jeans_mass_ref);
	
	total_mass_ref = give_ref (total_mass_max, total_mass_min);
	total_mass_err = give_err (total_mass_max, total_mass_ref);
	
	dm_dens_ref = give_ref (dm_dens_max, dm_dens_min);
	dm_dens_err = give_err (dm_dens_max, dm_dens_ref);
	
	bias_ref = give_ref (bias_max, bias_min);
	bias_err = give_err (bias_max, bias_ref);
	
	double chi_sq=0.;
	
	for (int i = 0; i < nBin; i++) {
	//	chi_sq+=pow(total_mass_ref.read(i)-jeans_mass_ref.read(i),2.)/
	//	(pow(total_mass_err.read(i),2.)+ pow(jeans_mass_err.read(i),2.));
		
	}
	
	printf("chi_sq = %g\n", chi_sq/nBin);
	
	
	
	
	ifstream inpFile;
	//inpFile.open("new_beta13_Man24_dm_dens_err.txt", ios::in);
	inpFile.open("new_beta13_Man24_dm_dens_err.txt", ios::in);
	//inpFile.open("new_betaZero_Man24_obs_gas_dens.dat", ios::in);
	array raio, dens, dens_err; raio.init(nBin-1); dens.init(nBin-1); dens_err.init(nBin-1);
	double rad, val, val_err;
	int i =0;
	while (inpFile >> rad>> val>> val_err) {
		if	(i!=2){
			raio.set(i, rad/oneMeterkpc);
			dens.set(i, val*(atom*1.e6));
			dens_err.set(i, val_err);
		}
	//	dm_dens_ref.set(i, val);//*(atom*1.e6));
	//	dm_dens_err.set(i, val_err);//*(atom*1.e6));
		i++;
	}
	inpFile.close();
	
	//ofstream outFile;
	//string dm_gamma =  "Mon30_gamma_13_err.txt";
	//outFile.open(dm_gamma.c_str(), ios::out);
	
	gas_mass_gauss = Gamma (raio, dens);
	
	writeArrayToFile(DoubleTimesArray(oneMeterkpc ,raio), gas_mass_gauss, "_dm_gamma_13");
	
	
	
	//ifstream inpFile;
	//inpFile.open("new_beta13_Man24_dm_dens_err.txt", ios::in);
	inpFile.open("new_betaZero_Man24_dm_dens_err.txt", ios::in);
	//inpFile.open("new_betaZero_Man24_obs_gas_dens.dat", ios::in);
	//array raio, dens, dens_err; raio.init(nBin-1); dens.init(nBin-1); dens_err.init(nBin-1);
	//double rad, val, val_err;
	i =0;
	while (inpFile >> rad>> val>> val_err) {
		if	(i!=2){
			raio.set(i, rad/oneMeterkpc);
			dens.set(i, val*(atom*1.e6));
			dens_err.set(i, val_err);
		}
		//	dm_dens_ref.set(i, val);//*(atom*1.e6));
		//	dm_dens_err.set(i, val_err);//*(atom*1.e6));
		i++;
	}
	inpFile.close();
	
	//ofstream outFile;
	//string dm_gamma =  "Mon30_gamma_13_err.txt";
	//outFile.open(dm_gamma.c_str(), ios::out);
	
	gas_mass_gauss = Gamma (raio, dens);
	
	writeArrayToFile(DoubleTimesArray(oneMeterkpc ,raio), gas_mass_gauss, "_dm_gamma_ZERO");
	
	
	
	
	
//	double * dm_dens_par; dm_dens_par = new double[8];
//	double guess_parameters[4] = {dm_dens_ref.read(2), 1.55373e22, 1., 3.};
//	dm_dens_par = fit_double_power_law_with_error_return_fit_and_error(radius_obs, 
//																	   dm_dens_ref, dm_dens_err, guess_parameters);
	//dm_dens_par = fit_density_PL_w_error_return_parameters_set_alpha_beta(radius_obs, dm_dens_ref, dm_dens_err, guess_parameters);
	
//	dm_dens_par = fit_density_PL_w_error_return_parameters_set_alpha_beta(raio, dens, dens_err, guess_parameters);
	
//	dm_dens_par = fit_density_PL_w_error_return_parameters_set_rs(radius_obs, dm_dens_ref, dm_dens_err, guess_parameters);

//	dm_dens_par = fit_double_power_law_with_error_return_fit_and_error(raio, dens, dens_err, guess_parameters);
																	   
	
//printf("(rho0, a, alpha, beta) = (%g, %g, %g, %g)  and (err1, err2, err3, err4) = (%g, %g, %g, %g) \n",
//		   dm_dens_par[0], dm_dens_par[1], dm_dens_par[2], dm_dens_par[3], dm_dens_par[4], dm_dens_par[5], dm_dens_par[6], dm_dens_par[7]);
	//dm_dens_ref = massFromDensityAssumeSpherical(radius_obs, dm_dens_ref);
	//dm_dens_err = massFromDensityAssumeSpherical(radius_obs, dm_dens_err);
	
//*/	
	
/*	ifstream inpFile, inpFile2;
	inpFile.open("new_betaZero_Man24_total_mass_err.txt", ios::in);
	inpFile2.open("new_betaZero_Man24_obs_gas_dens.dat", ios::in);
	double true_mass[nBin], true_mass_err[nBin];
	//double gas_mass[nBin], gas_mass_err[nBin];
	double rad, val, val_err;
	int i =0;
	while (inpFile >> rad>> val>> val_err) {
		true_mass[i]= val*oneMsunkg;
		true_mass_err[i] = val_err*oneMsunkg;
		i++;
	}
	inpFile.close();
	while (inpFile2 >> rad>> val>> val_err) {
		gas_mass_gauss.set(i,val);
		gas_mass_err.set(i, val_err);
		i++;
	}
	inpFile2.close();
	
	gas_mass_gauss = massFromDensityAssumeSpherical(radius_obs, gas_mass_gauss);
	gas_mass_err = massFromDensityAssumeSpherical(radius_obs, gas_mass_err);
	
	for (int i = 0; i < nBin; i++) {
		dm_dens_ref.set(i,true_mass[i]-gas_mass_gauss.read(i));
		//dm_dens_err.set(i,0.);
	}
	//dm_dens_err = massFromDensityAssumeSpherical(radius_obs, dm_dens_err);
	for (int i = 0; i < nBin; i++) {
		double value = pow(pow( true_mass_err[i]/ true_mass[i],2.) + pow( gas_mass_err.read(i)/gas_mass_gauss.read(i),2.),.5)*dm_dens_ref.read(i); 
		printf("val = %g A %g, B %g\n", value, true_mass_err[i], gas_mass_err.read(i));
		dm_dens_err.set(i, value);
	}
//*/	
	string name_he_err_2 = name_file+"_he_err.txt";
	outFile0.open(name_he_err_2.c_str(), ios::out);
	
	string name_jeans_err_2 = name_file+"_jeans_mass_err.txt";
	outFile1.open(name_jeans_err_2.c_str(), ios::out);
	
	string name_total_err_2 = name_file+"_GAS_MASS_err.txt";
	outFile2.open(name_total_err_2.c_str(), ios::out);
	
	string name_bias_err_2 = name_file+"_bias_err.txt";
	outFile3.open(name_bias_err_2.c_str(), ios::out);

	
	string name_dm_dens_err_2 = name_file+"_dm_MASS_err.txt";
	outFile4.open(name_dm_dens_err_2.c_str(), ios::out);
	
	
	
	for (int j = 0; j < nBin; j++) {
		outFile0 << radius_obs.read(j)*oneMeterkpc << " " << he_mass_ref.read(j)/oneMsunkg << " " << he_mass_err.read(j)/oneMsunkg<< endl;
		//outFile1 << radius_obs.read(j)*oneMeterkpc << " " << jeans_mass_ref.read(j)/(atom*1.e6) << " " << jeans_mass_err.read(j)/(atom*1.e6)<< endl;
		//outFile4 << radius_obs.read(j)*oneMeterkpc << " " << dm_dens_ref.read(j)/(atom*1.e6) << " " << dm_dens_err.read(j)/(atom*1.e6)<< endl;
		outFile4 << radius_obs.read(j)*oneMeterkpc << " " << dm_dens_ref.read(j)/oneMsunkg << " " << dm_dens_err.read(j)/oneMsunkg<< endl;
	
		outFile1 << radius_obs.read(j)*oneMeterkpc << " " <<  gas_mass_gauss.read(j)/oneMsunkg << " " << jeans_mass_err.read(j)/oneMsunkg<< endl;
		//outFile1 << radius_obs.read(j)*oneMeterkpc << " " <<  jeans_mass_ref.read(j)/oneMsunkg << " " << jeans_mass_err.read(j)/oneMsunkg<< endl;
		outFile2 << radius_obs.read(j)*oneMeterkpc << " " <<  total_mass_ref.read(j)/oneMsunkg << " " << total_mass_err.read(j)/oneMsunkg<< endl;
		//outFile2 << radius_obs.read(j)*oneMeterkpc << " " << total_mass_ref.read(j)/(atom*1.e6) << " " << total_mass_err.read(j)/(atom*1.e6)<< endl;
		
		outFile3 << radius_obs.read(j)*oneMeterkpc << " " << bias_ref.read(j) << " " << bias_err.read(j)<< endl;
		
	}
	outFile0.close();
	outFile1.close();
	outFile2.close();
	outFile3.close();
	outFile4.close();
	
	
	
	
	
	
	

//	total_mass_max_out = save_max_value(total_mass_max_out, total_mass_gauss);
//	if (nGauss==0) total_mass_min_out = CopyArray(total_mass_gauss);
//	else total_mass_min_out = save_min_value(total_mass_min_out, total_mass_gauss);
	
	
	
//	jeans_mass_max_out = save_max_value(jeans_mass_max_out, jeans_mass_gauss);
//	if (nGauss==0) jeans_mass_min_out = CopyArray(jeans_mass_gauss);
//	else jeans_mass_min_out = save_min_value(jeans_mass_min_out, jeans_mass_gauss);
	
	
/*	string name_he_err = name_file+"_total_err.txt";
	outFile0.open(name_he_err.c_str(), ios::out);
	string name_jeans_err = name_file+"_jeans_err.txt";
	outFile1.open(name_jeans_err.c_str(), ios::out);
	
	for (int j = 0; j < nBin; j++) {
		outFile0 << radius_obs.read(j)*oneMeterkpc << " " << total_mass_ref_out.read(j)/oneMsunkg << endl;
		//" " << total_mass_err_out.read(j)/oneMsunkg<< endl;
		outFile1 << radius_obs.read(j)*oneMeterkpc << " " << jeans_mass_ref_out.read(j)/oneMsunkg << endl;
		//" " << jeans_mass_err_out.read(j)/oneMsunkg<< endl;
	}
	outFile0.close();
	outFile1.close();
	
	
	
	
	
/*	jeans_mass_ref_out = give_ref (jeans_mass_max_out, jeans_mass_min_out);
	jeans_mass_err_out = give_err(jeans_mass_max_out, jeans_mass_ref_out);
	
	total_mass_ref_out = give_ref (total_mass_max_out, total_mass_min_out);
	total_mass_err_out = give_err(total_mass_max_out, total_mass_ref_out);
	
	
	string name_he_err = name_file+"_total_err.txt";
	outFile0.open(name_he_err.c_str(), ios::out);
	string name_jeans_err = name_file+"_jeans_err.txt";
	outFile1.open(name_jeans_err.c_str(), ios::out);
	*/
	
	return 0;
	
}
