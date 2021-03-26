/*
 *  clusterMakeItHappen.h
 *  
 *
 *  Created by Catarina Silva Fernandes  2010.
 *  Free of use! 
 *
 */

/* "Enough is enough!" Samuel Jackson (on "Snakes on a plane")
 */

#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <locale>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "clusterSetParameters.h"
#include "arrayMathematics.h"
#include "arrayFunctionsForCluster.h"

using namespace std;

double coutVariation(array Init, array Final) { //return a double maximum variation between correspondent members
	int i = Init.numberBins();
	
	double diference[2];
	diference[0] = 0.;
	
	for (int j=i*0.1; j < i*0.7; j++) {
		diference[1] = abs(( Init.read(j) - Final.read(j))/Final.read(j));
		
		if (abs(diference[1]) > abs(diference[0])) {
			diference[0] = diference[1];
		}
		
	}
	return diference[0];
} // return the maximum diference between Init and Final array


//When the parameters are defined/calculated nParameter == true
//until then nParameter == false. In this way we control functions that depend on others
void cluster::InitCalculusParameters () {
	nGasDensityInit = false;
	nGasGammaInit = false;
	nGasMassInit = false; 
	nGasDensityFinal = false;
	nGasGammaFinal = false;
	nGasMassFinal = false; 
	
	nDmDensityInit = false;
	nDmTildeInit = false;
	nDmGammaInit = false;
	nDmMassInit = false;
	nDmBetaInit = false; //
	nDmDensityFinal = false;
	nDmTildeFinal = false;
	nDmGammaFinal = false;
	nDmMassFinal = false;
	nDmBetaFinal = false; //
	
	nTotalMassInit = false;	
	nTotalMassHE = false;
	nTotalMassJeans = false;
	nTotalMassTrue = false;
	
	nDispVelInit = false;
	nDispVelFinal = false;
	
}
void cluster::InitParameters () { // all parameters set to zero
	
	printf("Starting all parameters, defining them as arrays, setting all to zero. Nbins = %d \n", nBin);	
	
	dm_disp_velocity = new array[1];dm_disp_velocity[0].init(nBin);
	dm_disp_velocity_error= new array[1]; dm_disp_velocity_error[0].init(nBin);
	dm_beta= new array[1]; dm_beta[0].init(nBin);
	dm_beta_error= new array[1]; dm_beta_error[0].init(nBin);
	dm_gamma= new array[1]; dm_gamma[0].init(nBin);
	dm_gamma_error= new array[1]; dm_gamma_error[0].init(nBin);
	Like2 = new array[1]; Like2[0].init(nBin);
	mass_nTh = new array[1]; mass_nTh[0].init(nBin);
	mass_HE = new array[1];	mass_HE[0].init(nBin);
	dm_mass = new array[1]; dm_mass[0].init(nBin);	
	dm_mass_error = new array[1]; dm_mass_error[0].init(nBin);
	dm_density = new array[1]; dm_density[0].init(nBin); 
	
	dm_density_error = new array[1];dm_density_error[0].init(nBin);
	mass_jeans = new array[1]; mass_jeans[0].init(nBin);
	mass_nTh = new array[1]; mass_nTh[0].init(nBin);
	gas_temperature = new array[1];	gas_temperature[0].init(nBin);
	gas_temperature_error= new array[1]; gas_temperature_error[0].init(nBin);
	gas_density = new array [1]; gas_density[0].init(nBin);
	gas_density_error = new array[1]; gas_density_error[0].init(nBin);
	gas_mass = new array[1]; gas_mass[0].init(nBin);
	gas_mass_error = new array[1]; gas_mass_error[0].init(nBin);
	mass_HE_error = new array[1]; mass_HE_error[0].init(nBin);
	
	
	YpCRs = new array[1]; YpCRs[0].init(nBin);
	MagneticDensity = new array[1]; MagneticDensity[0].init(nBin);
	MagneticField = new array[1]; MagneticField[0].init(nBin);
	
	GasDensityInit = new array[1]; GasDensityInit[0].init(nBin);
	GasGammaInit = new array[1]; GasGammaInit[0].init(nBin);
	GasMassInit = new array[1]; GasMassInit[0].init(nBin); 
	GasDensityFinal = new array[1]; GasDensityFinal[0].init(nBin);
	GasGammaFinal = new array[1]; GasGammaFinal[0].init(nBin);
	GasMassFinal = new array[1]; GasMassFinal[0].init(nBin);
	GasMassTesting = new array[1]; GasMassTesting[0].init(nBin); 
	GasDensityTesting = new array[1]; GasDensityTesting[0].init(nBin);
	GasGammaControl = new array[1]; GasGammaControl[0].init(nBin);
	
	DmDensityInit = new array[1]; DmDensityInit[0].init(nBin);
	DmTildeInit = new array[1]; DmTildeInit[0].init(nBin);
	DmGammaInit = new array[1]; DmGammaInit[0].init(nBin);
	DmMassInit = new array[1]; DmMassInit[0].init(nBin);
	DmBetaInit = new array[1]; DmBetaInit[0].init(nBin);
	DmDensityFinal = new array[1]; DmDensityFinal[0].init(nBin);
	DmTildeFinal = new array[1]; DmTildeFinal[0].init(nBin);
	DmGammaFinal = new array[1]; DmGammaFinal[0].init(nBin);
	DmMassFinal = new array[1]; DmMassFinal[0].init(nBin);
	DmBetaFinal = new array[1]; DmBetaFinal[0].init(nBin);
	DmBetaTesting = new array[1]; DmBetaTesting[0].init(nBin);
	DmTildeTesting = new array[1]; DmTildeTesting[0].init(nBin);
	DmDensityTesting = new array[1]; DmDensityTesting[0].init(nBin);
	DmMassTesting = new array[1]; DmMassTesting[0].init(nBin);
	
	TotalMassInit = new array[1]; TotalMassInit[0].init(nBin);
	TotalMassHE = new array[1]; TotalMassHE[0].init(nBin);
	TotalMassTrue = new array[1]; TotalMassTrue[0].init(nBin);
	TotalMassTrueTesting = new array [1]; TotalMassTrueTesting[0].init(nBin);
	TotalMassTesting = new array[1];  TotalMassTesting[0].init(nBin);
	TotalMassJeansInit = new array[1]; TotalMassJeansInit[0].init(nBin);
	TotalMassJeansTesting = new array[1]; TotalMassJeansTesting[0].init(nBin);
	TotalMassJeansFinal = new array[1]; TotalMassJeansFinal[0].init(nBin);
	
	TotalMassBiasCRs = new array[1]; TotalMassBiasCRs[0].init(nBin);
	TotalMassBiasB = new array[1];TotalMassBiasB[0].init(nBin);
	MassBiasCRsMax = new array[1];MassBiasCRsMax[0].init(nBin);
	MassBiasCRsMin = new array[1]; MassBiasCRsMin[0].init(nBin);
	MassBiasTurb1 = new array[1]; MassBiasTurb1[0].init(nBin);
	MassBiasTurb2 = new array[1]; MassBiasTurb2[0].init(nBin);
	MassBiasTurb3 = new array[1]; MassBiasTurb3[0].init(nBin);
	MassBiasBMin = new array [1]; MassBiasBMin[0].init(nBin);
	MassBiasBMax = new array[1]; MassBiasBMax[0].init(nBin);
	Temperature = new array[1];  Temperature[0].init(nBin);
	TemperatureTesting = new array[1];TemperatureTesting[0].init(nBin);
	TemperatureGamma = new array[1];TemperatureGamma[0].init(nBin);
	
	DispVelTesting = new array[1]; DispVelTesting[0].init(nBin);
	DispVelInit = new array[1]; DispVelInit[0].init(nBin);
	DispVelFinal = new array[1]; DispVelFinal[0].init(nBin);
	
		
	TotalDensityFinal = new array[1]; TotalDensityFinal[0].init(nBin);

	Cnonthermal = new array[nConfig];
	Cjeans = new array [nConfig];
	for (int i = 0; i < nConfig; i++) {
		Cnonthermal[i].init(nBin);
		Cjeans[i].init(nBin);
	}
	//Cparameter[0].init(nBin);

	CircularV = new array[1]; CircularV[0].init(nBin);
	
	// Pode ser que eu queira inserir outros profiles para a velocidade
	RadialAcceleration = new array[1]; RadialAcceleration[0].init(nBin);
	v_r = new array[1]; v_r[0].init(nBin);
	v_theta = new array[1]; v_theta[0].init(nBin);
	
}
void cluster::InitCluster () { // defines radial scale, init parameters and calculus
	
	//Initializes all parameters as an array and sets values in bin to Zero
	InitParameters();
	
	// SETS RADIUS, in LOG scale respectively linear numbers
	printf("Defining radius in Ln, Log10, Lin and /R200 = %g (pc) \n", 1.e3*R200);
	
	//parameters_fit_dm_mass = new array [1]; parameters_fit_dm_mass[0].init(8);
	
	RLn = new array[1]; RLn[0].init(nBin);
	RLog10 = new array[1]; RLog10[0].init(nBin);
	RLin = new array[1]; RLin[0].init(nBin);
	Rkpc = new array[1]; Rkpc[0].init(nBin);
	
	RLn[0] = setLnScale(rMin, rMax, nBin);
	RLog10[0] = setLog10scale(rMin, rMax, nBin);	//OK
	RLin[0] = expPowerArray(RLn[0]);	//OK
	Rkpc[0] = ArrayDividesDouble(RLin[0], 1.e3*onePcMtr*R200);
	
	
	RLin_error = new array[1]; RLin_error[0].init(nBin);
	// When the terms are calculated it becomes true
	InitCalculusParameters(); 
	
}

void cluster::InitCluster_Limits(double r_min_kpc, double r_max_kpc){
	//Initializes all parameters as an array and sets values in bin to Zero
	InitParameters();
	
	SetRMin(r_min_kpc*1.e3*onePcMtr,1.e2); //the 2nd value is irrelevant...
	SetRMax(r_max_kpc*1.e3*onePcMtr, 1.e2); //1st value in m
	
	// SETS RADIUS, in LOG scale respectively linear numbers
	printf("Defining radius in Ln, Log10, Lin and /R200 = %g (pc) \n", 1.e3*R200);
	
	parameters_fit_dm_mass = new array [1]; parameters_fit_dm_mass[0].init(8);
	
	RLn = new array[1]; RLn[0].init(nBin);
	RLog10 = new array[1]; RLog10[0].init(nBin);
	RLin = new array[1]; RLin[0].init(nBin);
	Rkpc = new array[1]; Rkpc[0].init(nBin);
	
	RLn[0] = setLnScale(rMin, rMax, nBin);
	RLog10[0] = setLog10scale(rMin, rMax, nBin);	//OK
	RLin[0] = expPowerArray(RLn[0]);	//OK
	Rkpc[0] = ArrayDividesDouble(RLin[0], 1.e3*onePcMtr*R200);
	
	// When the terms are calculated it becomes true
	InitCalculusParameters(); 
	
	
}

void cluster::ResetCluster() {
	
	delete []	GasDensityInit;
	delete []   GasGammaInit;
	delete []	GasMassInit;
	
	delete []	DmDensityInit ;
	delete []	DmGammaInit ;
	delete []	DmMassInit ;
	delete []	DmBetaInit ;
	delete []	DmTildeInit ;
	
	delete []	TotalMassInit ;
	delete []	TotalMassInit;	
	delete []	TotalMassHE;
	delete []	TotalMassJeansInit;
	delete []	TotalMassTrue;
	delete []	TotalMassTesting;
	
	delete []	MassBiasCRsMax;
	delete []	MassBiasCRsMin;
	delete []	MassBiasTurb1;
	delete []	MassBiasTurb2;
	delete []	MassBiasTurb3;
	delete []	MassBiasBMax;
	delete []	MassBiasBMin;
	
	delete []	DispVelTesting;
	delete []	DispVelInit;
	delete []	DispVelFinal;
	
	delete []	Temperature ;
	delete []	YpCRs;
	
	delete []	RLn; 
	delete []	RLog10;
	delete []	RLin; 
	delete []   Rkpc;
	
	GasDensityInit = NULL;
	GasGammaInit = NULL;
	GasMassInit = NULL;
	
	DmDensityInit = NULL;
	DmGammaInit = NULL;
	DmMassInit = NULL;
	DmBetaInit = NULL;
	DmTildeInit = NULL;
	
	TotalMassInit = NULL;
	TotalMassHE = NULL;
	TotalMassJeansInit = NULL;
	TotalMassTrue = NULL;
	TotalMassTesting = NULL;
	
	MassBiasCRsMax = NULL;
	MassBiasCRsMin= NULL;
	MassBiasTurb1= NULL;
	MassBiasTurb2= NULL;
	MassBiasTurb3= NULL;
	MassBiasBMax= NULL;
	MassBiasBMin = NULL;
	
	DispVelInit = NULL;
	DispVelFinal = NULL;
	DispVelTesting = NULL;
	
	Temperature = NULL;
	
	YpCRs = NULL;
	RLn = NULL; 
	RLog10 = NULL;
	RLin = NULL;
	Rkpc = NULL;
}

void cluster::SetTurbPressure(string bulkMotion, bool bulkRotation) {
	/*
	 
	 
	
	 for (int i = 0; i < nBin; i++) { // defined a profile for v_r
	 if (RLin[0].read(i) < 1.e1 * onePcMtr) {
	 value = -0.01*count2;
	 count2++;
	 }
	 if (RLin[0].read(i) > 5.e5 * onePcMtr) {
	 value -= 0.0001;
	 }
	 if ((RLin[0].read(i) > 1.e1 * onePcMtr)  && (RLin[0].read(i) < 5.e5 * onePcMtr) ) {
	 value =  -1.e4* pow(count1*1.e-3,2) ; 
	 count1++;
	 }
	 BulkMotion[0].set(i,value);
	 }
	
	*/
	
	
	double value = 0.;
	int count1, count2;
	count1 = 0.;
	count2 = 0.;
	
	if (bulkMotion == "case1") { //center static, ext inflow
		for (int i = 0; i < nBin; i++) { // defined a profile for v_r
			
			if ((Rkpc[0].read(i) < 0.0001) && (value <0.) ){
				value = 0.;
			}
			
			if ((Rkpc[0].read(i) > 0.001) && (value <0.) ){
				value = 0.;
				count2++;
			}
			if ((Rkpc[0].read(i) < 0.01) && ((value >0.) || (value == 0.)) ){
				value = 0.;
			}
			
			if ((Rkpc[0].read(i) >.2) && (value > -100.e3) ) {
				//value =  -5.e2* pow(count1*3.e-2,2) ; 
				value =  -3.e3* pow(count1*5.e-3,2) ;
				count1++;
			}
			if ((Rkpc[0].read(i) > 0.7) && (value < -100.e3)){
				value+= -1.e-4;
			}
			
			
			
			v_r[0].set(i,value);
			
		}
		writeArrayToFile(Rkpc[0],  ArrayDividesArray( v_r[0], CircularV[0]), "case1km-1");
		cout << "wrote case1 to file " << endl;
	}
	
	value = -100.;
	int count3 =0;
	count2 = 0;
	count1 = 0;
	if (bulkMotion == "case2") { //center inflow, ext inflow
		for (int i = 0; i < nBin; i++) { // defined a profile for v_r
			
			if ((Rkpc[0].read(i) < 0.0001) && (value <0.) ){
				value = -0.9e4;
			}
			
			if ((Rkpc[0].read(i) > 0.001) && (value <0.) ){
				value = -0.9e4 + pow(count2*1.e6,1./2.);
				count2++;
			}
			if ((Rkpc[0].read(i) < 0.1) && ((value >0.) || (value == 0.)) ){
				value = 0.;
			}
			
			if ((Rkpc[0].read(i) >.2) && (value > -100.e3) ) {
				value =  -5.e2* pow(count1*3.e-2,2) ; 
				//value =  -3.e3* pow(count1*5.e-3,2) ;
				count1++;
			}
			if ((Rkpc[0].read(i) > 0.7) && (value < -100.e3)){
				value+= -1.e-4;
			}
			
			
			v_r[0].set(i,value);
			
		}
		
		/*
		 for (int i = 0; i < nBin; i++) { // defined a profile for v_r
		 if (RLin[0].read(i) < 1.e1 * onePcMtr) {
		 value = -0.01*count2;
		 count2++;
		 }
		 if (RLin[0].read(i) > 5.e5 * onePcMtr) {
		 value -= 0.0001;
		 }
		 if ((RLin[0].read(i) > 1.e1 * onePcMtr)  && (RLin[0].read(i) < 5.e5 * onePcMtr) ) {
		 value =  -1.e4* pow(count1*1.e-3,2) ; 
		 count1++;
		 }
		 BulkMotion[0].set(i,value);
		 }
		 
		 */
		
		
		writeArrayToFile(Rkpc[0], ArrayDividesArray( v_r[0], CircularV[0]), "case2km-1");
		cout << "wrote case2 to file " << endl;
	}
	value = 100.;
	count2= 0;
	count1= 0;
	if (bulkMotion == "case3") { //center outflow, ext inflow
		for (int i = 0; i < nBin; i++) { // defined a profile for v_r
			
			if ((RLin[0].read(i) < R200*1.e3 *1.e-3* onePcMtr)&& (value > 0.) ){
				value = 1.e5 - pow(count2*5.e6,1./2.);
				count2++;
			}
			if ((RLin[0].read(i) < R200*1.e3 *1.e-2* onePcMtr) && ((value <0.) || (value == 0.)) ){
				value = 0.;
			}
			if ((RLin[0].read(i) > R200*1.e3 *1.e-1* onePcMtr)  && (value > -110.e3) ) {
				value =  -5.e2* pow(count1*5.e-3,2) ; 
				count1++;
			}
			if ((RLin[0].read(i) > R200*1.e3*1.e-1 * onePcMtr) && (value < -110.e3) ){
				value -= 0.01*count1;
			}
			v_r[0].set(i,value);
			
		}
		writeArrayToFile(Rkpc[0], DoubleTimesArray(1.e-3, v_r[0]), "case3kms-1");
		cout << "wrote case3 to file " << endl;
	}
	cout << v_r << endl;
	if (bulkRotation == false) {
		v_theta[0] = SetAllToZero(nBin);
	}
	
	count2 = 0;
	if (bulkRotation == true){
		for (int i = 0; i < nBin; i++) { // defined a profile for v_t^2
			
			if ((Rkpc[0].read(i) < 0.01) ){
				value = pow(0.01* CircularV[0].read(i),1);
				count2++;
			}
			if ((Rkpc[0].read(i) > 0.01) ){
				value = pow(0.2* CircularV[0].read(i),1);
				count2++;
			}
			v_theta[0].set(i,pow(value,2));
			
		}
		writeArrayToFile(Rkpc[0], ArrayDividesArray( v_theta[0], ArraySQ(CircularV[0])), "Rotationkm2-2");
		cout << "wrote rotation to file " << endl;
	}
	
	RadialAcceleration[0] = SetAllToZero(nBin);
	
	
}	


void cluster::SetClusterFromFile(string fileName){
	
	ifstream file; 
	file.open(fileName.c_str(),ios::in);
	double radius, density, temperature;
	
	int i =0;
	while (file>>radius>>density>>temperature) {
		i++;
	}
	cout << i << endl;
	SetNBin(i); //Defines number of bins for Cluster

	//InitCluster();
}


void cluster::SetGasFromFile(string fileName){
	
	ifstream file; 
	file.open(fileName.c_str(),ios::in);
	double radius, density, temperature;
	
	int j =0;
	while (file >> radius >>density >> temperature) {
		RLin[0].set(j,radius);	
		GasDensityFinal[0].set(j,density);
		Temperature[0].set(j,temperature);
		j++;
	}
	file.close();
	
	GasMassFinal[0] = massFromDensityAssumeSpherical(RLin[0], GasDensityFinal[0]);
	Rkpc[0] = ArrayDividesDouble(RLin[0], 1.e3*onePcMtr*R200);
	printf("Sets Gas DensityFinal and Temperature >>>> \n");
	
}
void cluster::SetDispVelDMbeta(string fileName){
	
	ifstream file; 
	file.open(fileName.c_str(),ios::in);
	double radius, dispVelocity, dmBeta;
	
	int j =0;
	while (file >> radius >>dispVelocity >> dmBeta) {
		RLin[0].set(j,radius);		
		DispVelFinal[0].set(j,dispVelocity);
		DmBetaInit[0].set(j,dmBeta);
		j++;
	}
	DispVelFinal[0] = dispVelFromTemperature(Temperature[0], DmBetaInit[0]);
	file.close();
	printf("Sets DispVelocity and Dm Beta \n");
	
}


void cluster::WriteAllTXT(string nameTemp) { // writes to fileNameG+String.txt
	
	//File with TIMING :)
	string time_string;
	time_t time_now = time(NULL);
	struct tm * timeinfo;
	timeinfo = localtime (&time_now);
	
	time_string = datetime_to_string(*timeinfo, "%H%M");
	
	
	const char *nameFile[20];
	
	string *test;
	test = new string[17];
	
	test[0] = time_string+"gasInit.txt";
	test[1] = time_string+"dmInit.txt";
	test[2] = time_string+"gasFinal.txt";
	test[3] = time_string+"dmFinal.txt";
	
	ofstream outGasInit; //density, gamma
	ofstream outGasFinal; //density, gamma, mass, temperature, temperatureGamma
	ofstream outDmInit; //density, gamma, mass, beta, dispersion velocity
	ofstream outDmFinal; //density, gamma, mass, beta, dispersion velocity
	
	outGasInit.open(test[0].c_str(),ios::out);	
	outDmInit.open(test[1].c_str(),ios::out);
	outGasFinal.open(test[2].c_str(),ios::out);
	outDmFinal.open(test[3].c_str(),ios::out);
	
	// we can put this on the first line
	//	outGasInit << "radius(/R200)" << " " << "gasDens(atom/cm3)" << " " << "gasGamma" << endl;
	//	outGasFinal << "radius(/R200)" << " " << "gasDens(atom/cm3)" << " " << "gasGamma" << " " << "mass(M_sun)" << " " << "temp(keV)" << " " << "tempGamma" << endl;
	for(int j = 0; j < nBin; j++) {	
		outGasInit << Rkpc[0].read(j);
		outGasInit << " " << (GasDensityInit[0].read(j)*1.e-6)/atom ;
		outGasInit << " " << GasGammaInit[0].read(j);
		outGasInit << " " << GasMassInit[0].read(j)/oneMsunkg<< endl;
		
		outGasFinal << Rkpc[0].read(j);
		outGasFinal << " " << (GasDensityFinal[0].read(j)*1.e-6)/atom ;
		outGasFinal << " " << GasGammaFinal[0].read(j);
		outGasFinal << " " << GasGammaControl[0].read(j);
		outGasFinal << " " << GasMassFinal[0].read(j)/oneMsunkg;
		outGasFinal << " " << Temperature[0].read(j);
		outGasFinal << " " << TemperatureGamma[0].read(j);
		outGasFinal << " " << ((GasDensityFinal[0].read(j)*1.e-6))*pow(RLin[0].read(j),2) <<endl;
		
	}
	outGasInit.close(); 
	outGasFinal.close();
	
	
	//	outDmInit << "radius(/R200)" << " " << "dmDens(atom/cm3)" << " " << "dmGamma" << " " << "dmMass(M_sun)" << " " << "sigma2_r(m2s-2)" <<endl;
	//	outDmFinal << "radius(/R200)" << " " << "dmDens(atom/cm3)" << " " << "dmGamma" << " " << "dmMass(M_sun)" << " " << "sigma2_r(m2s-2)" <<endl;
	for(int j = 0; j < nBin; j++) {	
		outDmInit << Rkpc[0].read(j);
		outDmInit << " " << (DmDensityInit[0].read(j)*1.e-6)/atom ;
		outDmInit << " " << DmGammaInit[0].read(j);
		outDmInit << " " << DmMassInit[0].read(j)/oneMsunkg;
		outDmInit << " " << pow(DispVelInit[0].read(j),2.)<<endl;
		
		outDmFinal << Rkpc[0].read(j);
		outDmFinal << " " << (DmDensityFinal[0].read(j)*1.e-6)/atom ;
		outDmFinal << " " << DmGammaFinal[0].read(j);
		outDmFinal << " " << DmMassFinal[0].read(j)/oneMsunkg;
		outDmFinal << " " << pow(DispVelFinal[0].read(j),2.)<<endl;
		
	}
	outDmInit.close(); //cout << "wrote to File Initial properties of the dm before HE+other = Jeans" << endl;
	outDmFinal.close();	//cout << "wrote to File Initial properties of the dm after HE+other = Jeans" << endl;
	

	writeArrayToFile(Rkpc[0], DoubleTimesArray( 1.e-6/atom, TotalDensityFinal[0]), "tDensity");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble( TotalMassInit[0],oneMsunkg), "MInit");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(TotalMassHE[0],oneMsunkg), "MHE");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(TotalMassJeansInit[0],oneMsunkg),"MJeans");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(TotalMassJeansFinal[0],oneMsunkg), "MJeansFinal");
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(TotalMassTrue[0],oneMsunkg), "MTrue.txt");

	writeArrayToFile(Rkpc[0], IntTimesArray(100,MassBiasCRsMax[0]), "BiasCRsMax");
	writeArrayToFile(Rkpc[0], IntTimesArray(100,MassBiasCRsMin[0]),"BiasCRsMin");
	writeArrayToFile(Rkpc[0], IntTimesArray(100,MassBiasTurb1[0]), "BiasTurb1");
	writeArrayToFile(Rkpc[0], IntTimesArray(100,MassBiasTurb2[0]), "BiasTurb2");
	writeArrayToFile(Rkpc[0], IntTimesArray(100,MassBiasTurb3[0]),"BiasTurb3");
	writeArrayToFile(Rkpc[0], IntTimesArray(100,MassBiasBMax[0]),"BiasBMax");
	writeArrayToFile(Rkpc[0], IntTimesArray(100,MassBiasBMin[0]), "BiasBMin");
	
}


/* FUNCTIONS 
 * 
 */

/* Total Mass calculus in perspective
 */
void cluster::TotalMassCalculus(string totalMass){
	
	if (totalMass == "Init" ) {
		TotalMassInit[0] = ArrayPlusArray(GasMassInit[0],DmMassInit[0]);
		nTotalMassInit = true;
		cout << "M_t (M_sun) from M_dm+M_g = " << TotalMassInit[0].read(nBin*0.9)/oneMsunkg << endl;
	}
	
	if (totalMass == "HE") {
		TotalMassHE[0] = MassFromHE(RLin[0], GasDensityFinal[0], Temperature[0]);
		nTotalMassHE = true;
		cout << "M_t (M_sun) from HE = " << TotalMassHE[0].read(nBin*0.90)/oneMsunkg << endl;
	}
	
	if (totalMass == "CRs") {
		//cout << "TotalMassTesting = M_HE + Delta(CRs) " << endl;
		if (nTotalMassHE == false) {
			cout << "TotalMassHE = f(r,T,rho_G) after HE=Jeans" << endl;
			
			TotalMassHE[0] = MassFromHE(RLin[0], GasDensityFinal[0], Temperature[0]);
			nTotalMassHE = true;
		}
		TotalMassTesting[0] = ArrayPlusArray(TotalMassHE[0], MassFromCRsPressure(RLin[0], Temperature[0], YpCRs[0], TotalMassHE[0]));
		cout << "M_t (M_sun) from HE + CRs = " << TotalMassTesting[0].read(nBin*0.90)/oneMsunkg << endl;
	}
	
	if (totalMass == "Turb") {
		if (nTotalMassHE == false) {
			cout << "TotalMassHE = f(r,T,rho_G) after HE=Jeans" << endl;
			
			TotalMassHE[0] = MassFromHE(RLin[0], GasDensityFinal[0], Temperature[0]);
			nTotalMassHE = true; 
		}
		array massTurb;
		massTurb.init(TotalMassHE[0].numberBins());
		massTurb = MassFromTurbPressure(RLin[0], v_r[0], v_theta[0]);
		writeArrayToFile(Rkpc[0], ArrayDividesDouble(massTurb, oneMsunkg), "massTurbulence");
		TotalMassTesting[0] = ArrayPlusArray(TotalMassHE[0], massTurb);
		cout << "M_t (M_sun) from HE + vel = " << TotalMassTesting[0].read(nBin*0.90)/oneMsunkg << endl;
		printf("Max Variation %e \n", coutVariation(massTurb,TotalMassHE[0]));
	}
	
	if (totalMass == "B") {
		if (nTotalMassHE == false) {
			cout << "TotalMassHE = f(r,T,rho_G) after HE=Jeans" << endl;
			
			TotalMassHE[0] = MassFromHE(RLin[0], GasDensityFinal[0], Temperature[0]);
			nTotalMassHE = true;
		}
		TotalMassTesting[0] = ArrayPlusArray(TotalMassHE[0], MassFromMagPressure(RLin[0], GasDensityFinal[0], MagneticField[0]));
		cout << "M_t (M_sun) from HE and B = " << TotalMassTesting[0].read(nBin*0.90)/oneMsunkg << endl;
	}
	
	if (totalMass == "Jeans") {
		cout << DmTildeTesting[0].read(10) << " " << DispVelTesting[0].read(10) << endl;
		TotalMassJeansTesting[0] = MassFromJeans(RLin[0], DmDensityTesting[0], ArraySQ(DispVelTesting[0]), DmBetaTesting[0]);
		cout << "M_t (M_sun) from Jeans = " << TotalMassJeansTesting[0].read(nBin*0.90)/oneMsunkg << endl;
	}
	
	if (totalMass != "Init" && totalMass != "HE" && totalMass != "CRs" && totalMass !="Turb" && totalMass != "B" && totalMass != "Jeans") {
		cout << "something is wrong... total mass "<< totalMass << "  was not calculated "<<endl;
		cout << "make sure you write Init, HE, CRs, Turb, B or Jeans " << endl;
	}
} 

/*	Initial Functions
 */
void cluster::GasDensityModel(){	
	
	if (gasModel == "hernquist") { //Hernquist
		GasDensityInit[0] = DensityHernquist(RLin[0], gasRho0, gasA);
		cout << "GAS DENSITY MODEL INIT, HERNQUIST <<<<<< " << endl;
		nGasDensityInit = true;
	}	// Hernquist
	
	
	if (gasModel == "beta") { //BetaModel
		GasDensityInit[0] = DensityBetaModel(RLin[0], gasRho0, gasA, betaModel);
		cout << "GAS DENSITY MODEL INIT, BETA with beta = " << betaModel << " <<<<<< " << endl;
		nGasDensityInit = true;
	}
	
	if (gasModel == "sersic") {
		GasDensityInit[0] = DensitySersicModel(RLin[0], gasRho0, gasA, sersicModel);
		cout << "GAS DENSITY MODEL INIT, SERSIC with nu = " << sersicModel << " <<<<<< " << endl;
		nGasDensityInit = true;
		
	}
	if (gasModel == "nfw") { //BetaModel
		GasDensityInit[0] = DensityNFW(RLin[0], gasRho0, gasA);
		cout << "GAS DENSITY MODEL INIT, NFW " << endl;
		nGasDensityInit = true;
	}
	if (gasModel != "hernquist" && gasModel != "beta" && gasModel != "sersic" && gasModel!= "nfw") {
		cout << "Missing other gas density profile models - RE-DO" << endl;
		nGasDensityInit = false;
	}
	
}	// OK
void cluster::DmDensityModel() {
	
	if (dmModel == "hernquist") { // Hernquist
		DmDensityInit[0] = DensityHernquist (RLin[0],dmRho0, dmA);
		cout << "DM DENSITY MODEL INIT, HERNQUIST <<<<<< " << endl;
		nDmDensityInit = true;
	}
	if (dmModel == "nfw") { //NFW
		DmDensityInit[0] = DensityNFW(RLin[0], dmRho0, dmA);
		cout << "DM DENSITY MODEL INIT, NFW <<<<<< " << endl;
		nDmDensityInit = true;
	}
	
	
	if (dmModel != "hernquist" && dmModel != "nfw") {
		cout << "Missing other dm density profile models - RE-DO" << endl;
		nDmDensityInit = false;
	}
}	// OK

void cluster::MassInitCalculusSpherical () { // gas and Dm init Mass asks for initial conditions if not set already
	
	while (nGasDensityInit == false) {
		string model;
		cout << "choose an Initial profile for gas write: [hernquist],[nfw],[beta], [sersic] " ;
		cin >> model;
		SetGasInitModel(model);
		GasDensityModel();
		
	}
	while (nDmDensityInit == false) {
		string model;
		cout << "choose an Initial profile for gdm write: [hernquist],[nfw] (make other profiles) " ;
		cin >> model;
		SetDmInitModel(model);
		DmDensityModel();
		
	}
	
	if (nGasMassInit == false) {
		GasMassInit[0] = massFromDensityAssumeSpherical(RLin[0], GasDensityInit[0]);
		nGasMassInit = true;
		cout << " GasMassInit has been calculated assuming spherical symmetry  " << endl;
		cout << "M_g (M_sun) = " << GasMassInit[0].read(nBin*0.9)/oneMsunkg << endl;
	}	
	
	if (nDmMassInit == false) {
		DmMassInit[0] = massFromDensityAssumeSpherical(RLin[0], DmDensityInit[0]);
		nDmMassInit = true;
		cout << " DmMassInit has been calculated assuming spherical symmetry  " << endl;
		cout << "M_dm (M_sun) = " << DmMassInit[0].read(nBin*0.9)/oneMsunkg << endl;
	}
	
	if(nTotalMassInit==false) {
		TotalMassCalculus("Init");
	}
	cout << " Gas, Dm Mass Init and Total = sum Init  " << endl;
}

void cluster::GammaInitCalculus (){  //gas and Dm
	while (nGasDensityInit == false) {
		string model;
		cout << "choose an Initial profile for gas write: [hernquist],[nfw],[beta], [sersic] " ;
		cin >> model;
		SetGasInitModel(model);
		GasDensityModel();
		
	}
	while (nDmDensityInit == false) {
		string model;
		cout << "choose an Initial profile for dm write: [hernquist],[nfw] (make other profiles) " ;
		cin >> model;
		SetDmInitModel(model);
		DmDensityModel();
		
	}
	
	
	if (nGasGammaInit == false) {	// gas gamma is calculated
		GasGammaInit[0] = Derivative(RLog10[0], Log10Array(GasDensityInit[0])); //dLogDens/dLogR
		nGasGammaInit = true;
	}
	if (nDmGammaInit == false) {	// dm gamma is calculated
		DmGammaInit[0] = Derivative(RLog10[0], Log10Array(DmDensityInit[0])); //dLogDens/dLogR
		nDmGammaInit = true;
		printf("Dm  Gamma: inner slope = %g and outer slope = %g \n", DmGammaInit[0].read(nBin*.1), DmGammaInit[0].read(nBin*.9));
	}
	
	cout << " Gas and Dm gamma Init  " << endl;
}
array BetaFollowsGammaLinear (array Gamma) {
	int g = Gamma.numberBins();
	
	array C;	
	C.init(g);
	
	C = DoubleTimesArray( (-1./6.), ArrayPlusInt(Gamma,1));
	
	return C;
}
void cluster::DmBetaInitCalculus () { //only Dm
	printf("DM BETA Set parameters \n");
	if (nDmBetaInit == false) {
		
		if (betaFollowsGamma == true) {
			if (nDmGammaInit == false) {
				GammaInitCalculus();
			}
			DmBetaInit[0] = BetaFollowsGammaLinear(DmGammaInit[0]);
			
			for ( int i = 0; i < nBin; i++) {	// attention if beta = 1
				
				if (DmBetaInit[0].read(i) >= 1.) {
					cout << "BETA >= 1  at Bin # " << i << endl;
					cin >> choice;
				}
			}
		}
		
		if (betaFollowsGamma == false) {
			cout << "We are assuming beta == 0 everywhere <<<<<<<< " << endl;
			DmBetaInit[0] = SetAllToZero(nBin);
		}
		
		nDmBetaInit = true;
		cout << " Beta Init for DM  >>>>>> "  << endl;
	}
}


/* During HE = Jeans
 */

void cluster::GetGasGammaFromDmProfile(){
	if (nDmGammaInit == false) {	// dm gamma is calculated
			DmGammaInit[0] = Derivative(RLog10[0], Log10Array(DmDensityInit[0])); //dLogDens/dLogR
			nDmGammaInit = true;
	}
	if (nDmBetaInit == false) {
		DmBetaInitCalculus();
	}
	GasGammaControl[0] = ArrayDividesArray(ArrayPlusArray(DmGammaInit[0],
														  ArrayPlusArray(DoubleTimesArray(2.,DmBetaInit[0]), 
																		 ArrayPlusArray(DoubleTimesArray(2./3.,ArrayTimesArray(DmBetaInit[0],
																															   Derivative(RLog10[0],Log10Array(ArraySQ(DispVelFinal[0]))))),
																						DoubleTimesArray(2./3.,Derivative(RLn[0],DmBetaInit[0])))))  ,
										   IntMinusArray(1,DoubleTimesArray(2./3.,DmBetaInit[0])));
	
	printf("Gamma Control for gas Density = %g, DmGammaInit = %g, DmBetaInit = %g \n", GasGammaControl[0].read(nBin*0.5), DmGammaInit[0].read(nBin*0.5), DmBetaInit[0].read(nBin*0.5));
}

/* During Jeans = HE+otherTerms
 */
void cluster::DmBetaTestingCalculus () { //only Dm
	
	if (betaFollowsGamma == true) {
		
		DmBetaTesting[0] = BetaFollowsGammaLinear(Derivative(RLog10[0], Log10Array(DmDensityTesting[0])));
		
		for ( int i = 0; i < nBin; i++) {	// attention if beta = 1
			
			if (DmBetaTesting[0].read(i) >= 1.) {
				cout << "BETA >= 1  at Bin # " << i << endl;
				cin >> choice;
			}
		}
	}
	
	else {
		cout << "We are assuming beta == 0 everywhere <<<<<<<< " << endl;
		DmBetaTesting[0] = SetAllToZero(nBin);
	}
}
void cluster::DmTildeTestingCalculus () { //only Dm
	
	
	DmTildeTesting[0] = YyFromdLogYydLogX (RLin[0], 
										   IntTimesArray(2, DmBetaTesting[0]), 
										   DmDensityTesting[0].read(0), 
										   DmDensityTesting[0]);
	
	//	nDmTildeInit = true;
	
	cout << " Dm rho tilde Testing >>>>>> "  << endl;
	
}

/* Rdm non-thermal
 */
double * random_nonThermal_mass(array Radius, array GasMassFinal, array HeMass, double * x){
	int r = Radius.numberBins();
	
	array NTmass;
	NTmass.init(r);
	
	double *rdm;
	rdm = new double[4];
	for (int i = 0; i < 4; i++) {
		rdm[i] = x[i];
	}
	int trial = 0;
	bool rdm_done = false;
	
	while (rdm_done == false) {
		
		rdm_done = true;
		
		NTmass = DensityPowerLaw(Radius, rdm[0], rdm[1], rdm[2], rdm[3]);
		for (int i = 0; i < r; i++) {
			if (HeMass.read(i)+NTmass.read(i) < GasMassFinal.read(i)){
				rdm_done = false; /*does not accept random number set as it fails our constraint */
			}
		}
		
		if(rdm_done == false){
			rdm[0] = rdm[0]*cat_rng_uniform_double_zero_one();
			rdm[1] = rdm[1]*cat_rng_uniform_double_zero_one();
			rdm[2] = rdm[2]*cat_rng_uniform_double_zero_one();
			rdm[3] = rdm[3]*cat_rng_uniform_double_zero_one();
		}
		trial++;
	}
	printf("found new configuration for NonThermalMass after %d trials \n", trial);
	printf("rho0 %g r_s %g alpha %g beta %g \n", rdm[0], rdm[1], rdm[2], rdm[3]);
	
	return rdm;
}

