/*
 *  clusterSetParameters.h
 *  
 *  Created by Catarina Silva Fernandes 2010.
 *  Free to be used!
 *
 */

/* "The imagination response is one of the principle mechanims in behavior motivation"
 * David Huron
 */

#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include "class.h"

using namespace std;

int cluster::ReadNumberBins (array A){
	return A.numberBins();
}

double cluster::ReadAscaleGas(){
	return gasA;
}


void cluster::SetRMin(double value, double valueCorrect) {
	rMin = value;
	rMinCorrect = valueCorrect;
}
void cluster::SetRMax(double value, double valueCorrect){
	rMax = value;
	rMaxCorrect = valueCorrect;
}
void cluster::SetR200(double value) {
	R200 = value;
}
void cluster::SetNBin(double value) {
	nBin = value;
}
void cluster::SetNConfig(double value){
	nConfig = value;
}
void cluster::SetFileName(string fileName) {
	FileName = fileName;
}

//Sets the model that we want the gas and dm to follow:
// model = hernquist, nfw, beta, sersic. If the names are wrong you can still change it in shell
void cluster::SetGasInitModel(string model){
	gasModel = model;
}
void cluster::SetDmInitModel(string model) {
	dmModel = model;
}
void cluster::SetGasRho0(double value){
	gasRho0 = value;
}
void cluster::SetDmRho0(double value){
	dmRho0 = value;
}
void cluster::SetAScale(double value){
	gasA = value;
	dmA = value;
}

//Sets the parameters from which this models are dependent
void cluster::SetBetaModel(double beta){
	betaModel = beta;
}
void cluster::SersicModel(double sersic){
	sersicModel = sersic;
}

//If beta=f(gamma) -> choose = true; if beta = 0 -> chose = false;
void cluster::SetBetaFollowGamma(bool choose){
	betaFollowsGamma = choose;
}

//Set the parameters
void cluster::SetCRsPressure(double Y0, double Psi){
	PsiCRs = Psi;
	Y0CRs = Y0;
}
void cluster::SetCRsProfileChange(bool choose){
	cosmicProfileChange = choose; 
	
}

void cluster::SetMagnProfileChange(bool choose){
	magProfileChange = choose;
}
void cluster::SetMagneticParameters(double B0, double BShape){
	magB0 = B0;
	magShape = BShape;
}


