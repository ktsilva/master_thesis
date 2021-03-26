/*
 *  MyClusterData.h
 *  
 *
 *  Created by Catarina Silva Fernandes on 28/12/2010.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */



#include <sstream>
#include <ctime>
#include <locale>
#include "libraries/clusterRdmWalk.h"
#include "MyProgramJumpMteCarlo.h"



void take_data_SI_to_files_A2052(string write_to_file){
	
	
	// Cluster A2052 - data in Mpc, in keV, in g/mˆ3, 1eV = 1.1604505e4K, 1g= 1.e-3kg
	string MyFile;
	
	double delta_keV = 1.e3* 1.1604505e4;
	double delta_Mpc = 1.e6*onePcMtr;
	double delta_dens = 1.e-3;
	int nBin = 8;
	array radius, radius_error, gas_density, gas_density_error, gas_temperature, gas_temperature_error;
	radius.init(nBin); gas_density.init(nBin); gas_temperature.init(nBin);
	radius_error.init(nBin); gas_density_error.init(nBin); gas_temperature_error.init(nBin);
	
	MyFile = "a2052.dat";
	ifstream file; file.open(MyFile.c_str(),ios::in);
	double parameter;
	int j,i; j=0; i=0;
	while (file >> parameter) {
		printf("Reading File i=%d \n",i);
		if (j < nBin*1) {radius.set(i,parameter*delta_Mpc); }	
		if (j >= nBin*1 && j < nBin*2) {radius_error.set(i,parameter*delta_Mpc);}
		if (j >= nBin*2 && j < nBin*3) {gas_temperature.set(i,parameter*delta_keV);}
		if (j >= nBin*3 && j < nBin*4) {gas_temperature_error.set(i,parameter*delta_keV);}
		if (j >= nBin*4 && j < nBin*5) {
			if (gas_temperature_error.read(i) > parameter){
				gas_temperature_error.set(i,parameter*delta_keV);
			}
		}
		if (j >= nBin*5 && j < nBin*6){gas_density.set(i,parameter*delta_dens);}
		if (j >= nBin*6 && j < nBin*7) {gas_density_error.set(i,parameter*delta_dens);}
		if (j >= nBin*7 && j < nBin*8) {
			if (gas_density_error.read(i) < parameter*delta_dens){
				gas_density_error.set(i,parameter*delta_dens);
			}
			//gas_density_error.set(i,gas_density_error.read(i)+gas_density.read(i));
		}
		
		i++;
		if (i == nBin) {i = 0;}
		j++;
	}
	file.close();

	ofstream data_file;
	MyFile = write_to_file + "_obs_gas_dens.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
			data_file << radius.read(k) << " " << gas_density.read(k) << " " << gas_density_error.read(k) << endl;
	}
	data_file.close();
	
	MyFile = write_to_file + "_obs_radius.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " << radius_error.read(k) << endl;
	}
	data_file.close();
	
	
	
	MyFile = write_to_file + "_obs_gas_temp.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " << gas_temperature.read(k) << " " << gas_temperature_error.read(k) << endl;
	}
	data_file.close();
	
	array gas_mass, mass_HE;
	gas_mass.init(nBin); mass_HE.init(nBin);
	
	mass_HE = MassFromHE( radius,  gas_density,  gas_temperature);
	gas_mass = massFromDensityAssumeSpherical(radius, gas_density);
	
	MyFile = write_to_file + "_obs_massHE.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " << mass_HE.read(k)/oneMsunkg << endl;
	}
	data_file.close();

	MyFile = write_to_file + "_obs_gas_mass.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " << gas_mass.read(k) << endl;
	}
	data_file.close();
	
	
	radius.delete_array(); radius_error.delete_array(); gas_density.delete_array(); gas_density_error.delete_array();
	gas_temperature.delete_array(); gas_temperature_error.delete_array(); gas_mass.delete_array(); mass_HE.delete_array();
	
}

void take_SI_fit_temp(string write_to_file, string read_density, string fit_dens, double * guess_dens,
					  string read_temperature, string fit_temp, double * guess_temp){
	
	int nBin = 0;
	string MyFile;
	
	double r, par, par_err;
	ifstream read_file; read_file.open(read_density.c_str(),ios::in);
	while (read_file >> r >> par >> par_err) {
		nBin++;
	}
	read_file.close();
	
	array radius, radius_error, gas_temperature, gas_temperature_error;
	radius.init(nBin); gas_temperature.init(nBin);
	radius_error.init(nBin); gas_temperature_error.init(nBin);
	
	read_file.open(read_temperature.c_str(),ios::in);
	int j = 0;
	while (read_file >> r >> par >> par_err) {
		radius.set(j, r);
		gas_temperature.set(j, par);
		gas_temperature_error.set(j, par_err);
		j++;
	}
	read_file.close();
	
	ofstream data_file;
	if (fit_temp == "double_power_law") {
		printf("Fitting Gas Temperature to double power law, writes parameters \n");
		double beta_temp = 3.;
		double temp_guess_dPL[4] = {guess_temp[0],guess_temp[1],guess_temp[2], beta_temp};//{1.e7, 1.e21, -.5, beta};
		double * fit_temp_dPL; fit_temp_dPL = new double [8];
		fit_temp_dPL = fit_double_power_law_fix_beta(radius, gas_temperature, gas_temperature_error, temp_guess_dPL);
		
		MyFile = write_to_file + "_temp_dPL_parameters.dat";
		data_file.open(MyFile.c_str(),ios::out);
		for (int k=0; k<8; k++) {
			data_file << fit_temp_dPL[k] << endl;
		}
		data_file.close();

	}
}


void take_data_SI_to_files_A1689(string write_to_file){
	
	
	// Cluster A1689 - data in pixel, in keV, in cmˆ-3, 1eV = 1.1604505e4K, 1g= 1.e-3kg
	string MyFile;
	
	double delta_keV = 1.e3* 1.1604505e4;
	double delta_pix = 185*1.e3*onePcMtr * (0.492/60.);
	double delta_ne = atom*1.17e3;

	int nBin=0;
	
	double par, par_errSup, par_errInf;
	MyFile = "all_temp_data.xcm"; // in keV
	ifstream file; file.open(MyFile.c_str(),ios::in);
	
	while (file >> par>> par_errSup >> par_errInf){
		nBin++;
		//printf("nBin = %d\n", nBin);
		
	}
	nBin--; // we are neglecting the last bin because its value is too high
	file.close(); 
	
	array gas_temperature, gas_temperature_error; gas_temperature.init(nBin); gas_temperature_error.init(nBin);
	array gas_density, gas_density_error; gas_density.init(nBin); gas_density_error.init(nBin);
	array radius, radius_error; radius.init(nBin); radius_error.init(nBin);

	//Gas Temperature
	int j=0;
	file.open(MyFile.c_str(),ios::in); 
	while (file >> par>> par_errSup >> par_errInf) {
		if (j != 0) {
		gas_temperature.set(j-1, par*delta_keV);
		double dif_sup = abs(par-par_errSup);
		double dif_inf = abs(par-par_errInf);
		gas_temperature_error.set(j-1, max(dif_sup, dif_inf)*delta_keV);
		}
		if (j == nBin-1) break;
		j++;
		
	}
	file.close();
	
	//Gas Density
	j=0;
	MyFile = "all_norm_data.xcm";  // in cmˆ3
	file.open(MyFile.c_str(),ios::in);
	
	while (file >> par>> par_errSup >> par_errInf) {
		if (j!= 0) {
			if (j == nBin) break;
			//printf("dens = %g errSup = %g  errInf = %g \n", par, par_errSup, par_errInf);
			
			gas_density.set(j-1, par*delta_ne);
			double dif_sup = par_errSup;//abs(par-par_errSup);
			double dif_inf = par_errInf;//abs(par-par_errInf);
			gas_density_error.set(j, abs(gas_density.read(j-1)-(max(dif_sup, dif_inf)*delta_ne)));
			//	gas_density_error_sup.set(j, par_errSup*delta_ne);
			//	gas_density_error_inf.set(j, par_errInf*delta_ne);
			printf("dens = %g err = %g \n", gas_density.read(j-1), gas_density_error.read(j-1));
			
		}
		j++;
	}
	file.close();
	
	j=0;
	MyFile = "annuli.dat"; 
	file.open(MyFile.c_str(),ios::in);
	while (file >> par) {
		if (j!=0) {
			if (j==nBin) break;
			radius.set(j-1, par*delta_pix);
			
		}
		j++;
	}
	file.close();
	
	ofstream data_file;
	
	MyFile = write_to_file + "_obs_radius.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " << radius_error.read(k) << endl;
	}
	data_file.close();
	
	
	MyFile = write_to_file + "_obs_gas_dens.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " << gas_density.read(k) << " " << abs(gas_density_error.read(k)) <<  endl;
	}
	data_file.close();
	
	MyFile = write_to_file + "_obs_gas_temp.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " << gas_temperature.read(k) << " " << abs(gas_temperature_error.read(k)) << endl;
	}
	data_file.close();
	
	//β = 0.8, ρ0 = 1.6×10−30 g cm^−3, and R = 1200 kpc.
	double beta = 0.8;
	double rho0 = 1.57e-30 * 1.e3; //kg m^−3
	double a = 500*1.e3 *onePcMtr;
	
	MyFile = write_to_file + "_obs_beta.dat";
	data_file.open(MyFile.c_str(),ios::out);
	
	for (int k=0; k<nBin; k++) {
		data_file << radius.read(k) << " " <<  rho0 * pow((1. + pow(radius.read(k)/a,2.)),-(3./2.)*beta) << " " <<  0.1*rho0 * pow((1. + pow(radius.read(k)/a,2.)),-(3./2.)*beta) << endl;
	}
	data_file.close();
	
	
	/////////////

	
	radius.delete_array(); radius_error.delete_array(); gas_density.delete_array(); gas_density_error.delete_array();
	gas_temperature.delete_array(); gas_temperature_error.delete_array(); 
	
}


void take_SI_fit_density(string write_to_file, string read_density, string fit_dens, double * guess_dens){

	
	int nBin = 0;
	string MyFile;
	
	double r, par, par_err;
	ifstream read_file; read_file.open(read_density.c_str(),ios::in);
	while (read_file >> r >> par >> par_err) {
		nBin++;
	}
	read_file.close();
	
	array radius, radius_error, gas_density, gas_density_error;
	radius.init(nBin); gas_density.init(nBin);
	radius_error.init(nBin); gas_density_error.init(nBin); 
	
	read_file.open(read_density.c_str(),ios::in);
	int j = 0;
	while (read_file >> r >> par >> par_err) {
		radius.set(j, r);
		gas_density.set(j, par);
		gas_density_error.set(j, par_err);
		j++;
	}
	read_file.close();
	
	ofstream data_file;
	if (fit_dens == "beta") {
		printf("Fitting Gas Density to Beta Model, writes parameters \n");
		double dens_guess_beta[3] = {guess_dens[0],guess_dens[1],guess_dens[2]};//{3.e-17, 1.e21, .5};
		double * fit_dens_beta; fit_dens_beta = new double [6];
		fit_dens_beta = fit_beta_profile_with_error_return_fit_and_error(radius, gas_density, gas_density_error, dens_guess_beta);
		
		MyFile = write_to_file + "_dens_beta_parameters.dat";
		data_file.open(MyFile.c_str(),ios::out);
		for (int k=0; k<6; k++) {
			data_file << fit_dens_beta[k] << endl;
		}
		data_file.close();
		
	}
	
	// for A1689: fit_dens_beta[2] = 0.8; fit_dens_beta[0] = 1.6e-30;//gcm-3; fit_dens_beta[1] = 1200 *1.e3*onePcMtr;
	
	if (fit_dens == "sersic") {
		printf("Fitting Gas Density to sersic Model, writes parameters \n");
		double dens_guess_sersic[3] = {guess_dens[0],guess_dens[1],guess_dens[2]}; //guess_dens[3]};//{3.e-17, 1.e21, .5};
		double * fit_dens_sersic; fit_dens_sersic = new double [6];
		fit_dens_sersic = fit_sersic_profile_with_error_return_fit_and_error(radius, gas_density, gas_density_error, dens_guess_sersic);
		
		MyFile = write_to_file + "_dens_sersic_parameters.dat";
		data_file.open(MyFile.c_str(),ios::out);
		for (int k=0; k<6; k++) {
			data_file << fit_dens_sersic[k] << endl;
		}
		data_file.close();
		
	}
	
	if (fit_dens == "double_power_law") {
		printf("Fitting Gas Density to Double Power Law Model, writes parameters \n");
		double beta_dens = 4.;
		double dens_guess_dPL[4] = {guess_dens[0],guess_dens[1],guess_dens[2], beta_dens};//{3.e-17, 1.e21, .5, beta};
		double * fit_dens_dPL; fit_dens_dPL = new double [8];
		fit_dens_dPL = fit_double_power_law_fix_beta(radius, gas_density, gas_density_error, dens_guess_dPL);
		
		MyFile = write_to_file + "_dens_dPL_parameters.dat";
		data_file.open(MyFile.c_str(),ios::out);
		for (int k=0; k<8; k++) {
			data_file << fit_dens_dPL[k] << endl;
		}
		data_file.close();		
	}
	

}