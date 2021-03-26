/*
 *	clusterGasNThDM.h
 *
 *  Created by Catarina Silva Fernandes 2010.
 *  Free to use!
 *
 */

/* As folhas caem no Outono
 */

#include <fstream>
#include <cmath>
#include <time.h>
#include <string>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>	
#include "clusterMakeItHappen.h"

using namespace std;

/*GasDensityFromDM: it's perfect
 *returns GasDensityFinal[0] profile and Temperture[0] from a initial dm density profile and a pseudo initial gas density
 *in the end it also calculates TotalMassHE[0], TotalMassJeansInit[0] and GasGammaFinal[0]
 */
void cluster::GasDensityFromDM(){ // finds properties of the gas with dm
	
	// Defines initial conditions for gas and dm, and total mass
	GasDensityModel(); //GasDensity[0] defined 
	DmDensityModel();	//DmDensity[0] defined
	GammaInitCalculus(); //calculates gamma for gas and dark matter density
	MassInitCalculusSpherical(); // M_t[0] = M_dm[0] + M_g[0]
	printf( "Initial profiles for gas and Dm \n");
	
	
	array TotalMassTesting, GasMassTesting, TemperatureTesting, DispVelTesting;
	TotalMassTesting.init(nBin);
	GasMassTesting.init(nBin);
	TemperatureTesting.init(nBin);
	DispVelTesting.init(nBin);
	
	
	while (nDmTildeInit == false) { 
		if (nDmBetaInit == false){
			DmBetaInitCalculus();
		} //DmBetaInit needs to be calculated first
		
		DmTildeInit[0] = YyFromdLogYydLogX (RLin[0], IntTimesArray(2, DmBetaInit[0]), 
											DmDensityInit[0].read(0), DmDensityInit[0]);
		
		nDmTildeInit = true;
		
		printf("Rho Tilde for DM Init >>>>> \n");
		
		
	} //Calculates DmTildeInit and DmBetaInit from DmDensityInit
	
	GasMassFinal[0] = CopyArray(GasMassInit[0]);
	printf("Testing Gas Mass and Density \n");
	
	double var=1.;
	
	while (var>1.e-6) { // infinite loop
		
		//TotalMassTesting = dmMassInit + NEWgasMass
		TotalMassTesting = ArrayPlusArray(GasMassFinal[0], DmMassInit[0]);
		
		//calculates a TemperatureTesting
		DispVelTesting =  dispVelFromTotalMass(RLin[0], DmTildeInit[0], TotalMassTesting);
		TemperatureTesting = temperatureFromDM(DispVelTesting, DmBetaInit[0]);
		printf("Testing Dispersion Velocity  and Temperature from Inital DmTilde and TotalMass \n");
		
		GasMassTesting = GasMass_DmProfileAndHE(RLin[0], TemperatureTesting, TotalMassTesting, GasDensityInit[0].read(0), DmMassInit[0].read(nBin-1));
		
		//compare - saves comparison in "var" and saves NEW Profile
		var = coutVariation(GasMassFinal[0], GasMassTesting);
		
		GasMassFinal[0] = CopyArray(GasMassTesting);
		GasDensityFinal[0] = DensityFromMassSpherical(RLin[0], GasMassFinal[0]);
		
	}
	
	//final value from loop we have a gasDensity
	nGasMassFinal = true;
	nGasDensityFinal = true;

	//let us define the final temperature Tg=Tdm
	Temperature[0] = CopyArray(TemperatureTesting);
	DispVelFinal[0] = CopyArray(DispVelTesting);
	
	nTemperature=true;
	
	if (nGasGammaFinal == false || nTemperatureGamma == false) {
		
		// Gas Gamma Final
		if (nGasDensityFinal == true) {
			GasGammaFinal[0] = Derivative(RLog10[0], Log10Array(GasDensityFinal[0]));
			nGasGammaFinal = true;
			printf("Gas Gamma Final >>>>>>>>>>\n");
		}
		else {
			printf("You can not calculate gasGammaFinal without knowing gasGamma from HE = Jeans");
		}
		
		// Temperature Final
		if(nTemperature == true) {
			TemperatureGamma[0] = Derivative(RLog10[0], Log10Array(Temperature[0]));
			nTemperatureGamma = true;
			printf("Temperature Slope  >>>>>>>>>>\n");
		}
		else {
			printf("you can not calculate temperatureGammaFinal without knowing Gas temperature from HE = Jeans ");
		}

	}
	
	
	//just checking
	if (betaFollowsGamma == true){
		writeArrayToFile(Rkpc[0], Derivative(RLog10[0],Log10Array(ArraySQ(DispVelTesting))), "DispVelInitGammaBeta");
	}
	if (betaFollowsGamma == false) {
		writeArrayToFile(Rkpc[0], Derivative(RLog10[0],Log10Array(ArraySQ(DispVelTesting))), "DispVelInitGammaBetaZero");
	}
	
	//calculates TotalMassHE[0] = f(gasDensity, temperature);
	TotalMassHE[0] = MassFromHE(RLin[0], GasDensityFinal[0], Temperature[0]);
	nTotalMassHE = true;
	
	TotalMassJeansInit[0] = MassFromJeans(RLin[0], DmDensityInit[0], ArraySQ(DispVelTesting),DmBetaInit[0]);
	
	printf("variation Testing/HE, Testing/Jeans = %f, %f \n ",
		   coutVariation(TotalMassHE[0], TotalMassTesting), coutVariation( TotalMassTesting, TotalMassJeansInit[0]));
	
	//TestingWithTheoryFromTeedy
	GetGasGammaFromDmProfile();
	writeArrayToFile(Rkpc[0], GasGammaControl[0], "GasGammaTheory");
	CircularV[0] = ArraySQRT(DoubleTimesArray(G,ArrayDividesArray(TotalMassHE[0], RLin[0])));
} 

//Adding NonThermal mass to MassHE see output
void cluster::TotalMassFromHEandCRs () {
	
	// steps from one extreme to the other
	int iMax = 30;
	double CRsMin = 0.1;
	double CRsMax = 0.5;
	double deltaCRs = (CRsMax - CRsMin)/(iMax*1.);
	
	int jMax = 10;
	double expCRsMax = 0.;
	double expCRsMin = -0.7;
	double deltaExpCRs = (expCRsMax - expCRsMin)/(jMax*1.);
	
	//we are going to Write the extreme values the one with less and more bias from HE	
	array MassBiasTesting, MassBiasMax, MassBiasMin, TotalMassTesting;
	//array TotalMassCRsMax, TotalMassCRsMin;
	
	MassBiasTesting.init(nBin);
	MassBiasMax.init(nBin);
	MassBiasMin.init(nBin);
	TotalMassTesting.init(nBin);
	
	MassBiasMin = SetAllToDouble(nBin,100.);
	//TotalMassCRsMax.init(nBin);
	//TotalMassCRsMin.init(nBin);
	//TotalMassCRsMin = SetAllToDouble(nBin, 1.e40*oneMsunkg); // we have to set the 1st min to a ridiculous large number
	
	cout << "Introducing CRs into TotalMassHE  >>>>>>>>>>" << endl;
	
	int testMAX =0;
	int testMIN =0;
	
	//for (int i = 0; i < iMax; i++) {
	//	for (int j = 0; j < jMax; j++) {
			
			//SetCRsPressure(CRsMin + (i*deltaCRs), expCRsMin + (j*deltaExpCRs));
			SetCRsPressure(0.1, -0.3);
			//printf("%g %g \n", CRsMin + (i*deltaCRs), expCRsMin + (j*deltaExpCRs));
			YpCRs[0] = DoubleTimesArray(Y0CRs,
										ArrayPowerDouble(ArrayDividesDouble(RLin[0], 
																			gasA),
														 PsiCRs)
									);
			//cout << "("<<CRsMin + (i*deltaCRs)<<","<<expCRsMin + (j*deltaExpCRs)<<") min" << endl;
			//TotalMassCalculus("CRs"); //saves to TotalMassTesting[0]
			TotalMassTesting = ArrayPlusArray(TotalMassHE[0], MassFromCRsPressure(RLin[0], Temperature[0], YpCRs[0], TotalMassHE[0]));
			MassBiasTesting = MassBiasCalculus(TotalMassHE[0], TotalMassTesting);
			
		/*	for (int k = 0; k< nBin;k++) {
				if (abs(MassBiasTesting.read(k)) > abs(MassBiasMax.read(k))){
					testMAX+=1;
				}
				if (abs(MassBiasTesting.read(k)) < abs(MassBiasMin.read(k))){
					testMIN+=1;
				}
			}
			
			if (testMAX > nBin*0.8){
				//TotalMassCRsMax = CopyArray(TotalMassTesting[0]);
				MassBiasCRsMax[0] = CopyArray(MassBiasTesting);
				MassBiasMax = CopyArray(MassBiasTesting);
				yp0max = Y0CRs;
				psimax = PsiCRs;
				cout << "("<<yp0max<<","<<psimax<<") max" << endl;
			}
			if (testMIN > nBin*0.8) {
				//TotalMassCRsMin = CopyArray(TotalMassTesting[0]);
				MassBiasCRsMin[0] = CopyArray(MassBiasTesting);
				MassBiasMin = CopyArray(MassBiasTesting);
				yp0min = Y0CRs;
				psimin = PsiCRs;
				
				cout << "("<<yp0min<<","<<psimin<<") min" << endl;
			}
			
			testMAX= 0;
			testMIN =0;		
		}
	}*/
	cout << " Has found max and minimum contribution of CRs " << endl;
	
	cout << "(yp0, psi) max&min  (" << yp0max << ","<< psimax<< ")  (" << yp0min << "," <<psimin << ")"<< endl; 
	writeArrayToFile(Rkpc[0], DoubleTimesArray(100.,MassBiasTesting), "raiosCosmicosMin");
	
	
	SetCRsPressure(0.5, -0.7);
	//printf("%g %g \n", CRsMin + (i*deltaCRs), expCRsMin + (j*deltaExpCRs));
	YpCRs[0] = DoubleTimesArray(Y0CRs,
								ArrayPowerDouble(ArrayDividesDouble(RLin[0], 
																	gasA),
												 PsiCRs)
								);
	//cout << "("<<CRsMin + (i*deltaCRs)<<","<<expCRsMin + (j*deltaExpCRs)<<") min" << endl;
	//TotalMassCalculus("CRs"); //saves to TotalMassTesting[0]
	TotalMassTesting = ArrayPlusArray(TotalMassHE[0], MassFromCRsPressure(RLin[0], Temperature[0], YpCRs[0], TotalMassHE[0]));
	MassBiasTesting = MassBiasCalculus(TotalMassHE[0], TotalMassTesting);
	
	/*	for (int k = 0; k< nBin;k++) {
	 if (abs(MassBiasTesting.read(k)) > abs(MassBiasMax.read(k))){
	 testMAX+=1;
	 }
	 if (abs(MassBiasTesting.read(k)) < abs(MassBiasMin.read(k))){
	 testMIN+=1;
	 }
	 }
	 
	 if (testMAX > nBin*0.8){
	 //TotalMassCRsMax = CopyArray(TotalMassTesting[0]);
	 MassBiasCRsMax[0] = CopyArray(MassBiasTesting);
	 MassBiasMax = CopyArray(MassBiasTesting);
	 yp0max = Y0CRs;
	 psimax = PsiCRs;
	 cout << "("<<yp0max<<","<<psimax<<") max" << endl;
	 }
	 if (testMIN > nBin*0.8) {
	 //TotalMassCRsMin = CopyArray(TotalMassTesting[0]);
	 MassBiasCRsMin[0] = CopyArray(MassBiasTesting);
	 MassBiasMin = CopyArray(MassBiasTesting);
	 yp0min = Y0CRs;
	 psimin = PsiCRs;
	 
	 cout << "("<<yp0min<<","<<psimin<<") min" << endl;
	 }
	 
	 testMAX= 0;
	 testMIN =0;		
	 }
	 }*/
	cout << " Has found max and minimum contribution of CRs " << endl;
	
	cout << "(yp0, psi) max&min  (" << yp0max << ","<< psimax<< ")  (" << yp0min << "," <<psimin << ")"<< endl; 
	writeArrayToFile(Rkpc[0], DoubleTimesArray(100.,MassBiasTesting), "raiosCosmicosMax");
	
	
	writeArrayToFile(Rkpc[0], DoubleTimesArray(100.,MassBiasMin), "raiosCosmicosMin");
	writeArrayToFile(Rkpc[0], DoubleTimesArray(100.,MassBiasMax), "raiosCosmicosMax");

	

} //it's perfect!


void cluster::TotalMassFromHEandTurb () { //finds TotalMass[3], Temperature[3], GasDensity[3] = GasDensity[1]
	
	char vString[3];
	
	string v;
	array MassBiasTurbTesting;
	MassBiasTurbTesting.init(nBin);

	SetTurbPressure("case1",false);
	TotalMassCalculus("Turb");
	MassBiasTurbTesting = MassBiasCalculus( TotalMassHE[0], TotalMassTesting[0]);
	writeArrayToFile(Rkpc[0],DoubleTimesArray(100.,MassBiasCalculus( TotalMassHE[0], TotalMassTesting[0])),"BiasT1" );
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(MassFromTurbPressure(RLin[0], v_r[0], v_theta[0]), oneMsunkg), "TurbulentMass");
	SetTurbPressure("case2",false);
	TotalMassCalculus("Turb");
	MassBiasTurbTesting = MassBiasCalculus( TotalMassHE[0], TotalMassTesting[0]);
	writeArrayToFile(Rkpc[0],DoubleTimesArray(100.,MassBiasCalculus( TotalMassHE[0], TotalMassTesting[0])),"BiasT2" );
	writeArrayToFile(Rkpc[0], Derivative(RLin[0], v_r[0]), "dv_r2");
	
	
} //it's ok -  we can always get different velocity profiles


void cluster::TotalMassFromHEandB()  { //finds TotalMass[4], Temperature[4], GasDensity[4] = GasDensity[1]
	
	
	
	
	MagneticDensity[0] = CopyArray(GasDensityFinal[0]); //GasDensity[1] comes from HE(th)
	magRho0 = GasDensityFinal[0].read(0);
	
	// steps from one extreme to the other
	int iMax = 5;
	double B0Min = 2.e-6*Gauss;
	double B0Max = 30.e-6*Gauss;
	double deltaB0 = (B0Max - B0Min)/(iMax*1.);
	
	int jMax = 5;
	double BShapeMin = 0.1;
	double BShapeMax = 0.9;
	double deltaBShape = (BShapeMax - BShapeMin)/(jMax*1.);
	
	//we are going to Write the extreme values the one with less and more bias from HE	
	array MassBiasTesting, MassBiasMax, MassBiasMin, massBias;
	
	massBias.init(nBin);
	MassBiasTesting.init(nBin);
	MassBiasMax.init(nBin);
	MassBiasMin.init(nBin);
	MassBiasMin = SetAllToDouble(nBin,100.);
	cout << "Introducing B into TotalMassHE  >>>>>>>>>>" << endl;
	
	double testMAX =0.;
	double testMIN =0.;
	double test =0.;
	
	//double GasDensity
	for (int i = 0; i< nBin; i++) {
		if (RLin[0].read(i) > 0.465570*gasA && RLin[0].read(i) < 0.465572*gasA  ) {
			
			printf("hello ");
		}
	}
	
	MagneticField[0] = DoubleTimesArray(4.8e-6*Gauss,
										ArrayPowerDouble(ArrayDividesDouble(GasDensityFinal[0],
																			gasRho0), 
														 0.6));
	
	TotalMassCalculus("B");
	MassBiasTesting = MassBiasCalculus(TotalMassHE[0], TotalMassTesting[0]);
		   
	massBias = CopyArray(MassBiasTesting);
		

	for (int i = 0; i < iMax; i++) {
		for (int j = 0; j < jMax; j++) {
			
			SetMagneticParameters(B0Min+(i*deltaB0), BShapeMin+(j*deltaBShape)); //minimum bfield in GC
			
			MagneticField[0] = DoubleTimesArray(B0Min+(i*deltaB0),
												ArrayPowerDouble(ArrayDividesDouble(GasDensityFinal[0],
																					gasRho0), 
																 BShapeMin+(j*deltaBShape)));
			
			TotalMassCalculus("B");
			MassBiasTesting = MassBiasCalculus(TotalMassHE[0], TotalMassTesting[0]);
			
			for (int k = 0; k< nBin;k++) {
				
				test += MassBiasTesting.read(k);
				testMAX+= MassBiasMax.read(k);
				testMIN+= MassBiasMin.read(k);
			}
			
			//	if (testMAX > nBin*0.8){
			
			if(abs(test) > abs(testMAX)){
				MassBiasBMax[0] = CopyArray(MassBiasTesting);
				MassBiasMax = CopyArray(MassBiasTesting);
				B0max = magB0;
				Shapemax= magShape;
				cout << "("<<B0max<<","<<Shapemax<<") max  " << testMAX << " " << test << " " << testMIN << endl;
			}
			//	if (testMIN > nBin*0.8) {
			if (abs(test) < abs(testMIN)){
				MassBiasBMin[0] = CopyArray(MassBiasTesting);
				MassBiasMin = CopyArray(MassBiasTesting);
				B0min = magB0;
				Shapemin = magShape;
				cout << "("<<B0min<<","<<Shapemin<<") min" << endl;
			}
			
			testMAX= 0;
			testMIN =0;		
			
			
		}
	}
	cout << " Has found max and minimum contribution of B " << endl;
	cout << "(B0, shape) max&min  (" << B0max << ","<< Shapemax<< ")  (" << B0min << "," <<Shapemin << ")"<< endl; 
	
	//Quick preview on outcome
		writeArrayToFile(ArrayDividesDouble(Rkpc[0], 1.), MagneticField[0], "MagneticField");
		writeArrayToFile(Rkpc[0], DoubleTimesArray(100., MassBiasMax), "BiasMaxFromBfield");
		writeArrayToFile(Rkpc[0],DoubleTimesArray(100.,massBias), "BiasFromBfield");
		writeArrayToFile(Rkpc[0],DoubleTimesArray(100.,MassBiasMin), "BiasMinFromBfield");
	
	SetMagneticParameters(B0max, Shapemax); //minimum bfield in GC
	
	MagneticField[0] = DoubleTimesArray(magB0,
										ArrayPowerDouble(ArrayDividesDouble(GasDensityFinal[0],
																			GasDensityFinal[0].read(nBin*0.1)), 
														 magShape));
	
	TotalMassCalculus("B");
	
	writeArrayToFile(Rkpc[0], DoubleTimesArray(1./oneMsunkg, TotalMassTesting[0]), "totalMassMagnetic");

	
} //it's perfect!



/*DmDensityFromMassHE:
 *retuns DmDensityFinal[0] profile and DispVelFinal[0]=f(Temperature,beta) - does not change
 *in the end it also calculates TotalMassJeansFinal[0] = TotalMassJeansInit[0] and DmGammaFinal[0]
 */
void cluster::DmDensityFromMassHE(){
	
	//Let totalMassTrue be equal to totalMassFromHE
	TotalMassTrue[0] = CopyArray(TotalMassHE[0]);
	TotalDensityFinal[0] = DensityFromMassSpherical(RLin[0], TotalMassTrue[0]);
	
	writeArrayToFile(Rkpc[0], TotalDensityFinal[0], "totalDensity");
	writeArrayToFile(Rkpc[0], Derivative(RLog10[0], Log10Array(TotalDensityFinal[0])), "totalGamma");
	printf("Getting True Mass and Density ");
	
	//M_dm= M_t - M_gas & DmBeta profile is defined to be the same as initially
	DmMassFinal[0] = ArrayMinusArray(TotalMassTrue[0], GasMassFinal[0]);
	DmBetaFinal[0] = CopyArray(DmBetaInit[0]);
	
	printf("Is DmMass[max] = X GasMass[max]? X = %f \n", DmMassFinal[0].read(nBin*0.99)/GasMassFinal[0].read(nBin*0.99));
	
		   
	DmDensityFinal[0] = DensityFromMassSpherical(RLin[0], DmMassFinal[0]);
	
	nDmMassFinal = true;
	nDmDensityFinal = true;
	if (nDmDensityFinal == true) {
		DmGammaFinal[0] = Derivative(RLog10[0], Log10Array(GasDensityFinal[0]));
	}
	else {
		printf("you can not calculate Dm gamma without density");
	}
	
		   
	DispVelFinal[0] = dispVelFromTemperature(Temperature[0], DmBetaInit[0]);
	writeArrayToFile(Rkpc[0], DispVelFinal[0], "DispVelFinalTemperature");
	TotalMassJeansFinal[0] = MassFromJeans(RLin[0], DmDensityFinal[0], ArraySQ(DispVelFinal[0]),DmBetaInit[0]);
	printf("How far is MJeans from MHE? %f %% \n", coutVariation(TotalMassTrue[0], TotalMassJeansFinal[0])*100.);
	printf("How far is MJeansFinal from MJeansInitial (expect it to be zero)? %f %% \n", coutVariation(DmMassInit[0], DmMassFinal[0]) );
	
	
		   /*	Comparing DispVel From Mass and Temperature
			* 
			* 	DmBetaTesting[0] = CopyArray(DmBetaInit[0]);
			DmDensityTesting[0] = CopyArray(DmDensityFinal[0]);
			DmTildeTestingCalculus();
			DispersionVelocityTestingCalculusFromTotalMass(DmTildeTesting[0], TotalMassInit[0]);
			writeArrayToFile(Rkpc[0], DispVelTesting[0], "DispVelFinalMass");
			*/	
		   
		
		   
		   
		   //TotalMassJeansFinal[0] = CorrectExtremesFunction(RLog10[0], MassFromJeans(RLin[0], DmDensityFinal[0], ArraySQ(DispVelFinal[0]), DmBetaFinal[0]));
		   		   
	//	DmDensityFinal[0] = CorrectExtremesFunction(RLog10[0],DensityFromMassSpherical(RLin[0], DmMassFinal[0]));
	/*	
	 double var = 1.;
	 
	 while (var > 0.001) {
	 DmDensityTesting[0] = DensityFromMassSpherical(RLin[0], DmMassTesting[0]);
	 
	 if (betaFollowsGamma == false) {
	 DmBetaTesting[0] = SetAllToZero(nBin);
	 }
	 
	 if( betaFollowsGamma == true) {
	 DmBetaTestingCalculus(); //needs DmDensityTesting
	 }
	 
	 DispersionVelocityTestingCalculusFromTemperature(Temperature[0], DmBetaTesting[0]);
	 cout << "variation on disp vel should be zero when beta = 0  " << coutVariation(DispVelFinal[0], DispVelTesting[0]) << endl;
	 
	 DmMassFinal[0] = DmMass_GasProfileAndJeans( RLin[0], ArraySQ(DispVelTesting[0]), DmBetaTesting[0], TotalMassTrue[0], DmDensityInit[0].read(0), GasMassFinal[0].read(nBin-1));
	 
	 
	 // cout << "MDm/Mgas max = " << DmMassFinal[0].read(nBin-1)/GasMassFinal[0].read(nBin-1) << endl;
	 //cin >> choice;
	 
	 
	 DmDensityTesting[0] = DensityFromMassSpherical(RLin[0], DmMassFinal[0]);
	 
	 DmTildeTestingCalculus();
	 
	 TotalMassJeansTesting[0] = MassFromJeans(RLin[0], DmTildeTesting[0], ArraySQ(DispVelTesting[0]));
	 
	 var = coutVariation(TotalMassJeansTesting[0], TotalMassTrue[0]); 
	 cout <<" variation Dm Final and Testing " << var << endl;
	 
	 DmTildeTestingCalculus();
	 TotalMassCalculus("Jeans");
	 cout << "compare Jeans with HE (should be zero) " << coutVariation(TotalMassHE[0], TotalMassJeansTesting[0]) << endl;
	 cout << "compare Jeans with HE + other terms " << coutVariation(TotalMassTrue[0], TotalMassJeansTesting[0]) << endl;
	 
	 
	 DmMassFinal[0] = DmMass_GasProfileAndJeans( RLin[0], ArraySQ(DispVelTesting[0]), DmBetaTesting[0], TotalMassTrue[0], DmDensityInit[0].read(0), GasMassFinal[0].read(nBin-1));
	 
	 
	 DmMassTesting[0] = CopyArray(DmMassFinal[0]);
	 DispVelFinal[0] = CopyArray(DispVelTesting[0]);
	 }
	 */
	
	
	
	
	
}



//Assume that GasDensityFinal[0] and Temperature[0] are our profiles from observations
// we can make it to read from FILE for real observations
// First let us gneerate RDM configuration AFTER let us make an educated guess
void cluster::DmDensityCparameterFromMassTrue() {
	//We generate a rdm configuration OR we can make an educated guess
	
	// ROUND #1
	SetCRsPressure(0.46, -0.7);
	YpCRs[0] = DoubleTimesArray(Y0CRs,
								ArrayPowerDouble(ArrayDividesDouble(RLin[0], 
																	gasA),
												 PsiCRs)
								);
	TotalMassCalculus("CRs"); 
	TotalMassTrueTesting[0]= CopyArray(TotalMassTesting[0]);
	writeArrayToFile(Rkpc[0], ArrayDividesDouble(TotalMassTesting[0], oneMsunkg), "totalMassCRs") ;
	
	cout<< "total mass HE and True " << endl;
	cout << TotalMassHE[0].read(9000) << endl;
	cout << TotalMassTrueTesting[0].read(9000) << endl;
	TotalDensityFinal[0] = DensityFromMassSpherical(RLin[0], TotalMassTesting[0]); //total density Final-change to Testing if necessary
	
	writeArrayToFile(Rkpc[0], TotalDensityFinal[0], "totalDensity");
	writeArrayToFile(Rkpc[0], Derivative(RLog10[0], Log10Array(TotalDensityFinal[0])), "totalGamma");
	cout << "Getting True Mass and Density "<<endl;
	
	
	DmMassTesting[0] = ArrayMinusArray(TotalMassTrueTesting[0], GasMassFinal[0]); //M_dm = T_total - M_gas
	DmBetaFinal[0] = CopyArray(DmBetaInit[0]); // I will not change beta profile with dm density again
	
	
	DmDensityTesting[0] = DensityFromMassSpherical(RLin[0], DmMassTesting[0]);
	
	
	nDmMassFinal = true;
	nDmDensityFinal = true;
	DmGammaFinal[0] = Derivative(RLog10[0], Log10Array(DmDensityFinal[0]));
	
	
	DispVelTesting[0] = dispVelFromTemperature(Temperature[0], DmBetaInit[0]);
	DispVelFinal[0] = CopyArray(DispVelTesting[0]);
	writeArrayToFile(Rkpc[0], DispVelFinal[0], "DispVelFinalTemperature");
	
	
	TotalMassJeansTesting[0] = MassFromJeans(RLin[0], DmDensityTesting[0], ArraySQ(DispVelTesting[0]),DmBetaFinal[0]);
	//TotalMassJeansFinal[0] = CorrectExtremesFunction(RLog10[0], MassFromJeans(RLin[0], DmDensityFinal[0], ArraySQ(DispVelFinal[0]), DmBetaFinal[0]));
	cout << "Mtrue = MHe? " << coutVariation(TotalMassHE[0],TotalMassTrueTesting[0]) << endl;
	cout << "Mjeans = Mtrue? " << coutVariation(TotalMassTrueTesting[0], TotalMassJeansTesting[0]) << endl;
	cout << "variation dm Initial and Final " << coutVariation(DmMassInit[0], DmMassFinal[0]) << endl;
	
	
	
	
}



/* Read from file in C++
 
 	ifstream file; 
	file.open("TrueMassSIunits.txt",ios::in);
	

	double x[n],y[n];
	double xx, yy;
	
	while (file >> xx >>yy ) {
		if ((g>(N-n)/2) && (g <N-(N-n)/2)){
			x[g-500] = xx;
			y[g-500] = yy;
		}
		
		g++;
	}
	
	file.close();
	
*/