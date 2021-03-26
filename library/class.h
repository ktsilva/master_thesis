  /*
 *	class.h
 *
 * ? set constants from gsl MSKSA (meters, kilograms, seconds, amperes)
 * ? creates 2 classes 
 *
 *  Created by Catarina Silva Fernandes 2010.
 *  Free to use!
 *
 */

/* A neve cai no Inverno
 */

#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <gsl/gsl_const_mksa.h>

using namespace std;


double const pi = 3.14159;

// Constants in SI units we can convert things afterwards to more Astro "friendly" units
//Physics
double const c = GSL_CONST_MKSA_SPEED_OF_LIGHT; // speed of light in vacuum, c.
double const mu0 = GSL_CONST_MKSA_VACUUM_PERMEABILITY; // permeability of free space, \mu_0. 
double const epsilon0 = GSL_CONST_MKSA_VACUUM_PERMITTIVITY; // permittivity of free space, \epsilon_0. 
double const h = GSL_CONST_MKSA_PLANCKS_CONSTANT_H; // Planck's constant, h.
double const hbar = GSL_CONST_MKSA_PLANCKS_CONSTANT_HBAR; // Planck's constant divided by 2\pi, \hbar.
//double const Na = GSL_CONST_NUM_AVOGADRO; // Avogadro's number, N_a.
double const kB = GSL_CONST_MKSA_BOLTZMANN; // Boltzmann constant, k_B.
double const R0 = GSL_CONST_MKSA_MOLAR_GAS; // molar gas constant, R_0.
double const Gauss = GSL_CONST_MKSA_GAUSS; // magnetic field of 1 Gauss. 

//Astronomy&Astrophysics
double const G = GSL_CONST_MKSA_GRAVITATIONAL_CONSTANT; // gravitational constant, G. 
double const onePcMtr = GSL_CONST_MKSA_PARSEC; // distance of 1 parsec, pc. 
double const oneMsunkg = GSL_CONST_MKSA_SOLAR_MASS; // mass of the Sun. 
double const mu = .61;

//Atomic&Nuclear Physics
double const e = GSL_CONST_MKSA_ELECTRON_CHARGE; // charge of the electron, e.
double const oneEvJoule = GSL_CONST_MKSA_ELECTRON_VOLT; // energy of 1 electron volt, eV.
double const atom = GSL_CONST_MKSA_UNIFIED_ATOMIC_MASS; // unified atomic mass, amu. (atom)
double const mE = GSL_CONST_MKSA_MASS_ELECTRON; // mass of the electron, m_e.
double const mP = GSL_CONST_MKSA_MASS_PROTON; // mass of the protron, m_p.

//My stuff
double const densUmax = 95.71e-28; //kg/m³
double const densUmin = 93.49e-28;
double const Kgp = kB/(mu*mP*G);


// *************** START WITH CLASSES ************************* //

class array {
	
private:
	
	int nBin;	// number of entries in the array
	double * bin;	// in each bin we have a value set to double
	string nameFolder;	// name of all files
	
public:
	void setNameFolder (string name); 
	string readNameFolder();
	
	void init(int numberBins);	// initilizes the array with nBin = numberBins
	void delete_array();
	int numberBins ();
	void set(int Nbin, double value);	// sets a value for the parameter in bin
	double read(int bin);	// reads the value of the parameter in entry bin
	
	
};

// we initialize the array;
// we can read how many entries has the array and the value of one of the entries
// we can set a value (double) to an entry
void array::init (int numberBins) { 
	nBin = numberBins;
	
	bin = new double[nBin];
	
	// when initializing each entry with value double it sets all to zero
	for(int i = 0; i < nBin; i++){ 
		bin[i] = 0.;
	}
} 
void array::delete_array(){
	delete[] bin;
}
int array::numberBins() {
	return nBin;
}
void array::set (int Nbin, double value) {
	bin[Nbin] = value;
}
double array::read (int Nbin) {
	return bin[Nbin];
}

void array::setNameFolder (string name) {
	nameFolder = name;
}
string array::readNameFolder (){
	return nameFolder;
}


// we define our class as a collection of arrays
class cluster {
	
private:
	
	int nBin;	// number of bins in each array 
	int nConfig; //number of configurations we want to generate for M_nt
	double rMax, rMin;
	double rMaxCorrect, rMinCorrect;
	double R200;
	
	
	string gasModel;	// we set a model and parameters for Gas Density
	
	double gasRho0;
	double gasA;
	double betaModel;
	double sersicModel;
	
	string dmModel;	// we set a model and parameters for Dm Density
	double dmRho0;
	double dmA;
	
	double magB0;
	double magShape;
	double magRho0;
	double B0max, B0min, Shapemax, Shapemin;
	
	double PsiCRs;
	double Y0CRs;
	double yp0max, yp0min, psimax, psimin; // in the end we will know which parameters define max&min
	
	bool betaFollowsGamma;	//if DmBeta follows a pre-assigned relation w/ DmGamma
	bool betaChange;	// if DmBeta changes w/ loop == true
	bool magProfileChange;
	bool cosmicProfileChange;
	bool turbProfileChange;
	
public:
	
	array * parameters_fit_dm_mass;
	
	string FileName;
	
	//Arrays that define our cluster and intermidiate ones that we need to test and run
	array * RLn; //ln scale
	array * RLog10;	//log10 scale
	array * RLin; // 10^log10 scale (log scale in linear form so that phys.s is ok)
	array * Rkpc; // R/R200
	
	array * RLin_error;
	array * dm_disp_velocity, * dm_disp_velocity_error;
	array * dm_beta, * dm_beta_error;
	array * dm_gamma, * dm_gamma_error;
	array * dm_mass, * dm_mass_error;
	array * dm_density, * dm_density_error;
	array * mass_HE, * mass_HE_error;
	array * mass_jeans;
	array * mass_nTh;
	array * gas_temperature, * gas_temperature_error;
	array * gas_density, * gas_density_error;
	array * gas_mass, * gas_mass_error;
	array * Like2;
	
	array * YpCRs;
	array * MagneticDensity;
	array * MagneticField;
	
	array * GasDensityInit;
	array * GasGammaInit;
	array * GasMassInit;
	
	array * GasDensityFinal;
	array * GasGammaFinal;
	array * GasMassFinal;
	
	array * GasMassTesting;
	array * GasDensityTesting;
	
	array * GasGammaControl; //from theory to compare with our results
	
	bool   nGasDensityInit;
	bool   nGasGammaInit;
	bool   nGasMassInit;
	bool   nGasDensityFinal;
	bool   nGasGammaFinal;
	bool   nGasMassFinal;
	
	
	array * DmDensityInit;
	array * DmTildeInit;
	array * DmGammaInit;
	array * DmMassInit;
	array * DmBetaInit;
	
	array * DmDensityFinal;
	array * DmTildeFinal;
	array * DmGammaFinal;
	array * DmMassFinal;
	array * DmBetaFinal; 
	
	array * DmBetaTesting;
	array * DmTildeTesting;
	array * DmDensityTesting;
	array * DmMassTesting;
	
	array * Cnonthermal;
	array * Cjeans;
	array * TotalDensityFinal;
	
	bool  nDmDensityInit;
	bool  nDmTildeInit;
	bool  nDmGammaInit;
	bool  nDmMassInit;
	bool  nDmBetaInit;
	bool  nDmDensityFinal;
	bool  nDmTildeFinal;
	bool  nDmGammaFinal;
	bool  nDmMassFinal;
	bool  nDmBetaFinal;
	
	array * TotalMassInit;	
	array * TotalMassHE;
	array * TotalMassTrue;
	array * TotalMassTrueTesting;
	array * TotalMassBiasCRs;
	array * TotalMassBiasTurb;
	array * TotalMassBiasB;
	
	array * TotalMassJeansInit;
	array * TotalMassJeansTesting;
	array * TotalMassJeansFinal; 	
	
	array * MassBiasCRsMax;
	array * MassBiasCRsMin;
	array * MassBiasTurb1;
	array * MassBiasTurb2;
	array * MassBiasTurb3;
	array * MassBiasBMax;
	array * MassBiasBMin;
	
	array * TotalMassTesting;
	
	array * ConfigurationNT;
	
	bool  nTotalMassInit;	
	bool  nTotalMassHE;
	bool  nTotalMassJeans;
	bool  nTotalMassTrue;
	
	array * CircularV;
	
	array  * DispVelTesting;
	array  * DispVelInit;
	array  * DispVelFinal;
	
	bool   nDispVelInit;
	bool   nDispVelFinal;
	
	array * TemperatureTesting;
	
	array  * Temperature;
	array  * TemperatureGamma;
	bool   	nTemperature;
	bool   	nTemperatureGamma;
	
	array * v_r;
	array * v_theta;
	array * v_psi;
	array * RadialAcceleration;
	
	// We have to write here all functions that we will need for cluster
	int ReadNumberBins(array A); 
	double ReadAscaleGas();
	void SetRMin(double value, double valueCorrect);
	void SetRMax(double value, double valueCorrect);	
	void SetNBin(double value);
	void SetNConfig(double value);
	
	void SetFileName(string fileName);
	
	void SetGasInitModel(string model);
	void SetDmInitModel(string model);
	void SetGasRho0(double value);
	void SetDmRho0(double value);
	void SetAScale(double value);
	void SetR200(double value);
	void SetBetaModel(double beta);
	void SersicModel(double sersic);	
	void SetBetaFollowGamma(bool choose);
	void SetCRsPressure(double Y0, double Psi);
	void SetCRsProfileChange(bool choose);
	void SetTurbPressure(string bulkMotion, bool bulkRotation);
	void SetMagneticParameters(double B0, double BShape);
	void SetMagnProfileChange(bool choose);
	
	void InitCalculusParameters();
	void InitParameters();
	void InitCluster();
	void InitCluster_Limits(double r_min_kpc, double r_max_kpc);
	void ResetCluster();
	
	void WriteAllTXT(string name);	
	void InitClusterSetGasFromFile(string MyFile);
	
	void GasDensityModel();
	void DmDensityModel();
	void MassInitCalculusSpherical();
	void TotalMassCalculus(string totalMass);	
	void GammaInitCalculus();
	void DmBetaInitCalculus ();
	void DmTildeInitCalculus ();
	void GammaFinalCalculus();
	
	void DispersionVelocityTestingCalculusFromTotalMass (array dmTildeInit, array totalMassTesting);
	void TemperatureTestingCalculus (array dispVelTesting, array dmBetaInit);
	
	void DmBetaTestingCalculus();
	void DmTildeTestingCalculus ();
	void DispersionVelocityTestingCalculusFromTemperature(array temperature, array dmBetaTesting);
	void GetGasGammaFromDmProfile();
	
	void GasDensityFromDM();
	void TotalMassFromHEandCRs ();
	void TotalMassFromHEandTurb ();
	void TotalMassFromHEandB ();
	void DmDensityFromMassHE();
	void DmDensityCparameterFromMassTrue();
	
	void SetClusterFromFile(string MyFile);
	void SetGasFromFile(string MyFile);
	void SetDispVelDMbeta(string MyFile);
	
	array * NewConfiguration(double * parameters, int nBin);
	void Testing();
//	void SimpleGreedy(double rho0DM, double aDM, double alphaDM, double betaDM, double eR, double eA, double eP, int nRejection, int nGreedy);
	void SimpleGreedy1(double rhoNew, double aNew, double alphaNew, double betaNew);
	void MonteCarlo(double rho0DM, double aDM, double alphaDM, double betaDM, double e, int nMonteCarlo, int nJump);
	void Jwalk(int nJwalk);
	
	void ChangeNth(double rhoNew, double aNew, double alphaNew, double betaNew);
	void CompareGasDM(string nameFile);
	
	void SimmulatedAnneling(double density0, double aScale, double alpha, double beta);
	
	void Fit_dm_mass(double * guess);
	void Init_dm_properties_beta_zero(double * guess);
	void Non_thermal_mass_init(double r_s);
	void Jeans_mass_from_total_mass_true(array TotalMass, array GasMass_Fit, double * guess);
	void MonteCarlo(double * guess, int nRej, int nIter, double * step, double * init_conditions);
	void Err_dm_mass_assume_density_powerLaw(int nBin);
	void Err_HE_mass_from_data(array Density, array Err_density, array Temperature, array Err_temperature, array Mass_from_HE);
	void Err_dm_beta_from_gamma();
	
};



