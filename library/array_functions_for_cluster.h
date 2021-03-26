/*
 *  arrayFunctionsForCluster.h
 *
 *  Created by Catarina Silva Fernandes 2010.
 *  Free to use!
 *
 */

/* "I was born not knowing and have had only a little time to change that here and there."
 * Richard P. Feynman
 */

#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <locale>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
using namespace std;

// DENSITY PROFILES
array DensityPowerLaw(array radius, double density0, double scale, double alpha, double beta) {
	int r = radius.numberBins();
	array Density;
	Density.init(r);
	
	Density = DoubleDividesArray( 
								 density0, ArrayTimesArray ( ArrayPowerDouble(ArrayDividesDouble(radius, scale),alpha),
															ArrayPowerDouble(ArrayPlusInt( ArrayDividesDouble(radius, scale) ,1),(beta-alpha))));
	
	
	//rho = rho0 / ((r/a)^alpha * (1+ (r/a))^(beta-alpha))

	return Density;
}
array DensityHernquist(array radius, double density0, double scale) { 
	
	int r = radius.numberBins();
	array Density;
	Density.init(r);
	
	/*Density = DoubleDividesArray( 
								 density0, ArrayTimesArray ( ArrayDividesDouble(radius, scale),
															ArrayPowerInt(ArrayPlusInt( ArrayDividesDouble(radius, scale) ,1),3)));
	*/
	Density = DensityPowerLaw(radius, density0, scale, 1., 4.);
	
	//rho = rho0 / ((r/a) * (1+ (r/a))^3)
	
	return Density;
}	// OK // alpha = 1, beta = 4
array DensityBetaModel (array radius, double density0, double scale, double beta) {
	
	int r = radius.numberBins();
	
	array Density;
	Density.init(r);
	double power = -(3./2.)*beta;
	
	Density = DoubleTimesArray( 
							   density0, ArrayPowerDouble( 
														  ArrayPlusInt( 
																	   ArraySQ(ArrayDividesDouble(radius, scale)),1), 
														  power));
	
	
	 //rho = rho0 * (1 + (r/a)²))^(-3*beta/2)
	
	return Density;
}	// OK
array DensityNFW (array radius, double density0, double scale) {
	int r = radius.numberBins();
	
	array Density;
	Density.init(r);
	
	/*Density = DoubleDividesArray( density0, 
								 ArrayTimesArray(ArrayDividesDouble( radius, scale),
												 ArraySQ( ArrayPlusInt( ArrayDividesDouble(radius, scale),1 )))
								 );
	*/
	Density = DensityPowerLaw(radius, density0, scale, 1., 3.);
	
	//rho = rho0 * (1/( (r/a) * (1+(r/a))²))
	
	
	return Density;
}	// OK //alpha = 1, beta = 3
array DensitySersicModel (array radius, double density0, double scale, double shape) {
	
	int r = radius.numberBins();
	array Density;
	Density.init(r);
	
	double v = shape;
	double Aa = scale * pow(2, 1/v);
	double Pp = 0.5 * (1 - 0.6097 * v + 0.05563 * pow(v,2));
	
	Density = DoubleTimesArray(density0,
							   ArrayTimesArray( ArrayPowerDouble( ArrayDividesDouble(radius, Aa),Pp),
											   expPowerArray(IntTimesArray(-1,
																		   ArrayPowerDouble(ArrayDividesDouble(radius, Aa) ,v)) )));
	
	//rho = rho0 *(r/a')^p' * exp(-(r/a')^v)
	//p' = p/2
	//p = 1-0.6097v+0.05563v^2
	//a'= a*2^(1/v)
	
	return Density;
}

// MASS CALCULUS
/*	Mass from density profiles - spherical profile -  M(r) = Int_0^r 4*pi * rho * r'² dr'
 */
array massFromDensityAssumeSpherical (array radius, array density) {
	
	int r = radius.numberBins();
	array mass;
	mass.init(r);
	
	mass = IntegralTrapezeZeroToPointP (radius,
										DoubleTimesArray(4.*pi, ArrayTimesArray (density,ArraySQ(radius))));
	
	return mass;
	
}
array MassFromHE(array radius, array densityGas, array temperatureGas) {
	
	int r = radius.numberBins();
	
	array MHE;
	MHE.init(r);
	
	array DlogRho, DlogTemp;
	DlogRho.init(r);
	DlogTemp.init(r);
	
	DlogRho = Derivative(Log10Array(radius),
						 Log10Array(densityGas));
	DlogTemp = Derivative(Log10Array(radius),
						  Log10Array(temperatureGas));
	
	MHE = DoubleTimesArray( -Kgp, 
						   ArrayTimesArray(ArrayTimesArray(radius,
														   temperatureGas),
										   ArrayPlusArray(DlogRho, DlogTemp))
						   );
	DlogRho.delete_array(); DlogTemp.delete_array();
	return MHE;
}

array MassFromCRsPressure(array radius, array temperatureGas, array ratioPressureGasCRs, array massHE) {
	
	int r = radius.numberBins();
	
	array MCRs;
	MCRs.init(r);
	
	
	MCRs = ArrayPlusArray ( ArrayTimesArray( ratioPressureGasCRs,
											massHE), 
						   DoubleTimesArray ( -Kgp,
											 ArrayTimesArray(ArrayTimesArray(ArraySQ(radius),
																			 temperatureGas), 
															 Derivative(radius,
																		ratioPressureGasCRs)))
						   );
	return MCRs;
}
array MassFromMagPressure(array radius, array gasDensity, array MagneticField) {
	int r = radius.numberBins();
	
	array MB;
	MB.init(r);
	
	
	MB =  DoubleTimesArray ( -1./(2.*mu0*G), ArrayTimesArray( ArrayDividesArray(ArraySQ(radius),gasDensity), 
															 Derivative(radius, ArraySQ(MagneticField))));
	
	return MB;
	
}
array MassFromTurbPressure(array radius, array bulkMotion, array bulkRotation) {
	int r = radius.numberBins();
	
	array MTurb;
	MTurb.init(r);
	
	
	//MTurb = ArrayTimesArray( ArrayTimesArray( ArraySQ(radius),bulkMotion), 
	//								ArrayDividesDouble(Derivative(radius,bulkMotion),G));
	//	MTurb = ArrayTimesArray(ArrayDividesDouble(radius, G), ArraySQ(bulkMotion)  ); 
//	MTurb = ArrayMinusArray( ArrayTimesArray( ArrayTimesArray( ArraySQ(radius),
//															  bulkMotion), 
//											 //ArrayDividesDouble(Derivative(radius,bulkMotion),G)) ,
//											 ArrayDividesArray(bulkMotion,ArrayDividesDouble(radius,G))),
//							ArrayDividesDouble(ArrayTimesArray(radius,bulkRotation),G));	

//	MTurb = ArrayTimesArray(ArrayTimesArray(ArrayDividesDouble(ArraySQ(radius),G),
//											bulkMotion),
//							Derivative(radius,bulkMotion));
	MTurb = ArrayPlusArray( DoubleTimesArray( 1./G,ArrayTimesArray(radius, ArraySQ(bulkRotation))),
						   ArrayTimesArray(ArrayDividesDouble(ArraySQ(bulkMotion),G),
							radius));
											
	return MTurb;
}

array MassFromJeans(array radius, array dmDensity, array sigma2, array dmBeta) {
	int r = radius.numberBins();
	
	array M;
	M.init(r);
	M = ArrayDividesDouble( ArrayTimesArray( radius, 
											ArrayTimesArray( sigma2, ArrayPlusArray( IntTimesArray(2,dmBeta), ArrayPlusArray( Derivative(Log10Array(radius), Log10Array(dmDensity)) ,
																															 Derivative(Log10Array(radius), Log10Array(sigma2)))))),-G);
	return M;
}

//Density From Spherical 
array DensityFromMassSpherical(array radius, array Mass) {
	
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	C = DoubleTimesArray(1./(4.*pi), ArrayDividesArray(Derivative(radius, Mass),
													   ArraySQ(radius)));
	
	return C;						
	
}
array DensityFromMassSimple(array radius, array Mass){
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	for (int i = 0; i < r; i++) {
		double value;
		if (i==0) value = Mass.read(i) /((4./3.)*pi * pow(radius.read(i),3.));
		if (i !=0) {
			value = (3./(4.*pi))*abs(Mass.read(i)-Mass.read(i-1))/(pow(radius.read(i),3.)-pow(radius.read(i-1),3.));
		}
		C.set(i, value);
	}
	
	
	return C;
}

//Dm Density from DmTilde density and beta profile -1/6 (1+gamma_dm) very specific
array DmDensity_DmTildeBeta (array radius, array dmTilde) {
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	C = ArrayTimesArray(ArrayPowerDouble( ArrayDividesDouble( dmTilde,pow(dmTilde.read(0)*radius.read(0),1./3.)),3./2.),
						ArraySQRT(radius));
	
	return C;
}

// MASS
array MassBiasCalculus(array massEstimate, array massTrue) {
	int m = massTrue.numberBins();
	
	array C;
	C.init(m);
	
	C = ArrayDividesArray (ArrayMinusArray( massEstimate, massTrue), massTrue);
	
	return C;
}

//MassGas = 0.1 MasDm
array GasMassFromAlpha (array radius, array Alpha, double gasDens0Rdm, double dmMassMax) { //OK
	
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	array densityTest, massTest, massNormalised;
	densityTest.init(r);
	massTest.init(r);
	massNormalised.init(r);
	
	densityTest = YyFromdLogYydLogX (radius, Alpha, gasDens0Rdm, SetAllToInt(r,1)); 
	
	densityTest.set(0, densityTest.read(1));
	densityTest.set(r-1, densityTest.read(r-2));
	
	massTest = massFromDensityAssumeSpherical(radius, densityTest);
	
	double massTestMax = massTest.read(r-1);
	printf("CORRECTING GAS MASS TOWARDS Gas Mass MAX = 0.1 DM Mass MAX ");
	
	massNormalised = DoubleTimesArray( 0.1 *dmMassMax/massTestMax , massTest);
	
	C = CopyArray(massNormalised);
	
	densityTest.delete_array(); massTest.delete_array(); massNormalised.delete_array();
	return C;
}

array GasMass_DmProfileAndHE (array radius, array temperature, array totalMass, double rhoGasMin, double MassDmMax) {
	
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	
	array AlphaHE;
	
	AlphaHE.init(r);
	
	AlphaHE = ArrayPlusArray( IntTimesArray(-1,
											Derivative(Log10Array(radius), Log10Array(temperature)
													   )
											), 
							 IntTimesArray(-1,
										   DoubleTimesArray((G*mu*mP)/(kB),
															ArrayDividesArray(totalMass, 
																			  ArrayTimesArray(radius, 
																							  temperature
																							  )
																			  )
															)
										   )
							 );
	
	
	double gasDens0Rdm = rdmNumber(.1,10.)*rhoGasMin;
	
	C = GasMassFromAlpha(radius, AlphaHE, gasDens0Rdm, MassDmMax);
	//array GasDensityFromAlpha (array radius, array Alpha, double gasDens0Rdm, double dmDensMax)
	
	cout << " NEW GAS Mass >>>>>> " << C.read(r*0.9)/oneMsunkg << endl;
	AlphaHE.delete_array();
	
	return C;
}

//MassDm = 10 MassGas
array DmMassFromAlpha(array radius, array Alpha, double dmDens0Rdm, double gasMassMax) {
	
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	array densityTest, massTest, massNormalised;
	densityTest.init(r);
	massTest.init(r);
	massNormalised.init(r);
	
	densityTest = YyFromdLogYydLogX (radius, Alpha, dmDens0Rdm, SetAllToInt(r,1)); 
	
	densityTest.set(0, densityTest.read(1));
	densityTest.set(r-1, densityTest.read(r-2));
	
	massTest = massFromDensityAssumeSpherical(radius, densityTest);
	
	double massTestMax = massTest.read(r-1);
	printf("CORRECTING DM MASS TOWARDS GAS Mass MAX = 10 Gas Mass MAX ");
	
	massNormalised = DoubleTimesArray( 10.* gasMassMax/massTestMax, massTest);
	
	C = CopyArray(massNormalised);
	
	densityTest.delete_array(); massTest.delete_array(); massNormalised.delete_array();
	
	return C;
	
	
}

array DmMass_GasProfileAndJeans( array radius, array sigma2, array beta, array totalMass, double rhoDmMin, double MassGasMax) {
	
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	
	array AlphaJeans;
	
	AlphaJeans.init(r);
	
	AlphaJeans = ArrayPlusArray(IntTimesArray(-2, beta),
								ArrayPlusArray( IntTimesArray(-1,
															  Derivative(Log10Array(radius), Log10Array(sigma2)
																		 )
															  ), 
											   IntTimesArray(-1,
															 DoubleTimesArray((G),
																			  ArrayDividesArray(totalMass, 
																								ArrayTimesArray(radius, 
																												sigma2
																												)
																								)
																			  )
															 )
											   ));
	/*If The outcome will be a fucntion of DmDensity and Beta
	 AlphaJeans = ArrayPlusArray( IntTimesArray(-1,
	 Derivative(Log10Array(radius), Log10Array(sigma2)
	 )
	 ), 
	 IntTimesArray(-1,
	 DoubleTimesArray((G),
	 ArrayDividesArray(totalMass, 
	 ArrayTimesArray(radius, 
	 sigma2
	 )
	 )
	 )
	 )
	 );*/
	
	double dmDens0Rdm = rdmNumber(.1,10.)*rhoDmMin;
	
	C = DmMassFromAlpha(radius, AlphaJeans, dmDens0Rdm, MassGasMax);
	
	cout << " NEW DM Mass >>>>>> " << endl;
	AlphaJeans.delete_array();
	return C;
}

// Temperature from Vdisp and Beta Model
array temperatureFromDM (array dispVelTesting, array dmBeta){
	int d = dispVelTesting.numberBins();
	
	array C;
	C.init(d);
	
	C = DoubleTimesArray((mu*mP)/(kB), ArrayTimesArray(ArraySQ(dispVelTesting),IntMinusArray(1,DoubleTimesArray(2./3., dmBeta))));
	
	printf(" TEMPERATURE FROM DM >>>>>> \n");
	
	
	return C;
	
}


// Dispersion velocity
array dispVelFromTotalMass (array radius, array dmTildeInit, array totalMassTesting){
	int r = radius.numberBins();
	
	array C;
	C.init(r);
	
	array yAxis, integral;
	yAxis.init(r);
	integral.init(r);
	
	yAxis = ArrayDividesArray(ArrayTimesArray(dmTildeInit, totalMassTesting), ArraySQ(radius));
	integral = IntegralTrapezePointPtoInfty(radius, yAxis);
	
	C = ArraySQRT(ArrayTimesArray(DoubleDividesArray(G, dmTildeInit), 
												  integral));
	
	printf(" DISPERSION VELOCITY from Total Mass >>>>>> \n");
	
	yAxis.delete_array(); integral.delete_array();
	
	return C;
	
}

array dispVelFromTemperature(array temperature, array dmBetaTesting){
	int t = temperature.numberBins();
	
	array C;
	C.init(t);
	
	C = ArraySQRT(ArrayDividesArray( DoubleTimesArray( kB/(mu*mP),temperature), IntMinusArray(1, DoubleTimesArray(2./3., dmBetaTesting))));
	
	printf(" DISPERSION VELOCITY from Temperature >>>>>> \n");
	
	return C;
}


//correct Extremes
array CorrectExtremesFunction(array radiusLog10, array function){
	int f = function.numberBins();
	
	array C;
	C.init(f);
	C= CopyArray(function);
	double value;
	double deltaX = radiusLog10.read(1) - radiusLog10.read(0);
	array derivative;
	derivative.init(f);
	derivative=SetAllToZero(f);
	
	double derL, derR, der;
	derivative = Derivative(radiusLog10,Log10Array(function));
	for (int i = f*0.1; i >-1; i--){
		der = derivative.read(i);
		derL = derivative.read(i-2);
		derR= derivative.read(i+2);
		
		if (abs((derL-derR)/derR)> 1.e-2 || (i==0) || (i==1)){
			//value = C.read(i+1) - der*(radiusLog10.read(i+1)-radiusLog10.read(i-1))- derL*(radiusLog10.read(i-1)-radiusLog10.read(i-2));
			value = C.read(i+2) - (2.*deltaX*derR);
			//	value = (1/(radiusLog10.read(i+2)-radiusLog10.read(i+1))) *
			//			(C.read(i+1)*(radiusLog10.read(i+2)-radiusLog10.read(i+1)) +
			//				(C.read(i+2)-C.read(i+1))*(radiusLog10.read(i+1)-radiusLog10.read(i)) -
			//					2.*derR*(radiusLog10.read(i+1)+radiusLog10.read(i))*(radiusLog10.read(i+2)-radiusLog10.read(i+1)));
		}	
		C.set(i, value);
		
	}
	// ponto zero nao tem derivada à esquerda...
	for (int i = f*0.9;i<f;i++) {
		der = derivative.read(i);
		derL = derivative.read(i-2);
		derR= derivative.read(i+2);
		
		if (abs((derL-derR)/derL)> 1.e-2 || (i != f-1) || (i!= f-2)){
			//value = C.read(i-2) + der*(radiusLog10.read(i+1)-radiusLog10.read(i-1)) + derR*(radiusLog10.read(i+2)-radiusLog10.read(i+1));
			value = (2.*deltaX*derL) +C.read(i-2);
		}
		C.set(i, value);
	}
	
	derivative.delete_array();
	
	return C;
}


void writeParametersToFile (int nBins, double Rho0[], double aScale[], double alpha[], double beta[], string nameFile){
	
		//write with TIME
	string time_string;
	time_t time_now = time(NULL);
	struct tm * timeinfo;
	timeinfo = localtime (&time_now);
	time_string = datetime_to_string(*timeinfo, "%H%M");
	string myFile;
	
	myFile = time_string + nameFile + ".txt";
	
	ofstream outFile;
	outFile.open(myFile.c_str(),ios::out);
	
	for (int i = 0; i < nBins; i++) {
		
		outFile << i << " " << Rho0[i] << " " << aScale[i] << " " << alpha[i] << " " << beta[i] << endl;
		
	}
	
	outFile.close();
	
}


