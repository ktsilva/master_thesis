/*
 *  arrayErrorPropagation.h
 *  
 *
 *  Created by Catarina Silva Fernandes on 2010.
 * 
 *
 */

#include <time.h>
#include <fstream>
#include <sstream>
#include <cmath>
#include <ctime>
#include <locale>
#include <iostream>
#include <stdio.h>
#include <string.h>


array error_propagation_mass_HE( array mass_HE, array radius, array temp, array temp_err, array dens, array dens_err){

	int nBin = mass_HE.numberBins();
	
	array error; error.init(nBin);
	array d_rho, d_T; d_rho.init(nBin); d_T.init(nBin);
	d_rho = Gamma(radius, dens); d_T = Gamma(radius, temp);
	
	array dM_dT, dM_dgamma_dens, dM_dgamma_temp; dM_dT.init(nBin); dM_dgamma_dens.init(nBin); dM_dT.init(nBin);
	array section_rho, section_T; section_rho.init(nBin); section_T.init(nBin);
	
	dM_dT = Derivative(temp, mass_HE);
	dM_dgamma_dens = Derivative(d_rho, mass_HE);
	dM_dgamma_temp = Derivative(d_T,mass_HE);
	
	for (int i = 0 ; i < nBin; i++) {
		double value_rho, value_T;
		if (i == 0) {
			value_rho = pow(dens_err.read(i)/dens.read(i),2.) + 
			((pow(dens_err.read(i),2.)+pow(dens_err.read(i+1),2.))/pow(dens.read(i+1)-dens.read(i),2.));
			
			value_T = pow(temp_err.read(i)/temp.read(i),2.) + 
			((pow(temp_err.read(i),2.)+pow(temp_err.read(i+1),2.))/pow(temp.read(i+1)-temp.read(i),2.));
		}
		if (i == nBin-1) {
			value_rho = pow(dens_err.read(i)/dens.read(i),2.) + 
			((pow(dens_err.read(i),2.)+pow(dens_err.read(i-1),2.))/pow(dens.read(i)-dens.read(i-1),2.));
			
			value_T = pow(temp_err.read(i)/temp.read(i),2.) + 
			((pow(temp_err.read(i),2.)+pow(temp_err.read(i-1),2.))/pow(temp.read(i)-temp.read(i-1),2.));
		}
		if (i != 0 && i!= nBin-1) {
			value_rho = pow(dens_err.read(i)/dens.read(i),2.) + 
			((pow(dens_err.read(i+1),2.)+pow(dens_err.read(i-1),2.))/pow(dens.read(i+1)-dens.read(i-1),2.));
			
			value_T = pow(temp_err.read(i)/temp.read(i),2.) + 
			((pow(temp_err.read(i+1),2.)+pow(temp_err.read(i-1),2.))/pow(temp.read(i+1)-temp.read(i-1),2.));
		}
		
		section_rho.set(i, value_rho);
		section_T.set(i, value_T);
	}
	
	
	error = ArraySQRT(ArrayPlusArray ( ArraySQ(ArrayTimesArray(dM_dT, temp_err)), 
							ArrayPlusArray( ArrayTimesArray(ArraySQ(ArrayTimesArray( dM_dgamma_dens,d_rho)), section_rho) , 
										   ArrayTimesArray(ArraySQ(ArrayTimesArray( dM_dgamma_temp,d_T)), section_T))));
	coutOneArrayToScreen(ArrayDividesArray(error, mass_HE));
	printf("Temp\n");
	coutOneArrayToScreen(ArraySQ(ArrayTimesArray(dM_dT, temp_err)));
	printf("gammaDens\n");
	coutOneArrayToScreen(ArrayTimesArray(ArraySQ(ArrayTimesArray( dM_dgamma_dens,d_rho)), section_rho));
	printf("gammaTemp\n");
	coutOneArrayToScreen(ArrayTimesArray(ArraySQ(ArrayTimesArray( dM_dgamma_temp,d_T)), section_T));
	

	
	coutThreeArrayToScreen (mass_HE, error, ArrayDividesArray(error, mass_HE));
	
	return error;
	
}




array error_propagation_beta_model(array beta_function, array radius, double * fit){
	
	int nBin = beta_function.numberBins();
	array error; error.init(nBin);
	double rho0, r_c, beta, rho0_error, r_c_error, beta_error;
	rho0 = fit[0]; r_c = fit[1]; beta = fit[2]; rho0_error = fit[3]; r_c_error = fit[4]; beta_error = fit[5];
	for (int i =0; i < nBin; i++) {
		
		double value = pow(pow(beta_function.read(i)* rho0_error/rho0,2.) + 
						   pow(3*beta*beta_function.read(i)*(pow(radius.read(i)/r_c,2.)) *
							   (r_c_error/ (r_c*(pow(r_c,2.) + pow(radius.read(i),2.)))),2.)+
						   pow(beta*beta_function.read(i) * abs(log(1.+ pow(radius.read(i)/r_c,2.))*(-3./2.))*
							   beta_error,2.),.5);
		error.set(i, value);
	}
	
	return error;
	
}

array error_propagation_dPL(array doublePL_function, array radius, double *fit){
	
	int nBin = doublePL_function.numberBins();
	array error; error.init(nBin);
	double rho0, r_s, alpha, beta, rho0_error, r_s_error, alpha_error, beta_error;
	rho0 = fit[0]; r_s = fit[1]; alpha = fit[2], beta = fit[3];
	rho0_error = fit[4]; r_s_error = fit[5]; alpha_error = fit[6]; beta_error = fit[7];
	
	
	for (int i =0; i < nBin; i++) {
		
		double value = pow(pow(doublePL_function.read(i)* rho0_error/rho0,2.) + 
						   pow(doublePL_function.read(i)* r_s_error *
							   (alpha + (beta*radius.read(i)/r_s))/(r_s + radius.read(i)) ,2.)+
						   pow( doublePL_function.read(i) * alpha_error *
							   abs(log((radius.read(i)+r_s)/radius.read(i))),2.)+
						   pow( doublePL_function.read(i) * beta_error *
							   abs(log(r_s/(radius.read(i)+r_s))),2.),.5);
		error.set(i, value);
	}
	
	return error;
	
}

array error_propagation_one_var(array function, array var1, array var1_error) {
	
	int nBin = function.numberBins();
	
	array error; error.init(nBin);
	
	error = ArraySQRT( ArrayTimesArray (ArraySQ(Derivative(var1, function)),
										ArraySQ(var1_error)));
	return error;
}

/*array error_propagation_mass_HE(array function, array var1, array var1_error) {
	
	int nBin = function.numberBins();
	
	array error; error.init(nBin);
	
	error = ArrayTimesArray( ArrayDividesArray(function, var1),var1_error);
	
	return error;
}

*/

array error_propagation_two_var(array function, array var1, array var1_error, array var2, array var2_error) {
	
	int nBin = function.numberBins();
	
	array error; error.init(nBin);
	
	error = ArraySQRT(ArrayPlusArray( ArrayTimesArray (ArraySQ(Derivative(var1, function)),
													   ArraySQ(var1_error)) ,
									 ArrayTimesArray (ArraySQ(Derivative(var2, function)),
													  ArraySQ(var2_error))));
	
	
	
	for (int i = 0; i < nBin; i++) {
		printf("%d (%g, %g) // ((%g*%g, %g*%g))\n",i, ArrayTimesArray (ArraySQ(Derivative(var1, function)),
												ArraySQ(var1_error)).read(i), ArrayTimesArray (ArraySQ(Derivative(var2, function)),
								ArraySQ(var2_error)).read(i), Derivative(var1, function).read(i), var1_error.read(i),
			   Derivative(var2, function).read(i), var2_error.read(i));
	}
	return error;
}



array err_density_powerLaw (double * fit_parameters, double * r, int nBin) {

	// Model power-law = rho0 /( pow(r/r_s,alpha) * pow(1+r/r_s,beta-alpha)) 
	
	double rho0 = fit_parameters[0]; double r_s =  fit_parameters[1]; double alpha =  fit_parameters[2]; double beta =  fit_parameters[3];
	double err_rho0 =  fit_parameters[4]; double err_r_s =  fit_parameters[5]; double err_alpha =  fit_parameters[6]; double err_beta =  fit_parameters[7];
	
	double density[nBin];
	for (int i = 0; i < nBin; i++) density[i] = rho0 * pow((r_s + r[i])/r[i],alpha) * pow(r_s/(r_s+r[i]),beta); 
	
	double j0[nBin],j1[nBin], j2[nBin], j3[nBin];
	for (int i = 0; i < nBin; i++) {
		j0[i] = density[i]/rho0;
		j1[i] = density[i] * (1./(r_s + r[i]))* (alpha + ((beta * r[i])/r_s));
		j2[i] = density[i] * log((r[i]+r_s)/r[i]);
		j3[i] = density[i] * log(r_s/(r[i]+r_s));
	}
	array error;
	error.init(nBin);
	
	for (int i = 0 ; i < nBin ; i++) error.set(i,pow( pow(j0[i] * err_rho0,2.) + pow (j1[i] * err_r_s,2.) 
								+ pow(j2[i]*err_alpha,2.)+ pow(j3[i]*err_beta,2.),0.5));
	
	return error;
}

void cluster::Err_dm_mass_assume_density_powerLaw(int nBin) {
	
	// Model mass from density as power-law = Integral[4.*pi*rho0*pow(r_s,beta)*pow(r[i],(2.-alpha))*pow((r[i]+r_s),(alpha-beta))]

	double rho0 = parameters_fit_dm_mass[0].read(0);double r_s =  parameters_fit_dm_mass[0].read(1);
	double alpha = parameters_fit_dm_mass[0].read(2);double beta =  parameters_fit_dm_mass[0].read(3);

	double err_rho0 = parameters_fit_dm_mass[0].read(4);double err_r_s = parameters_fit_dm_mass[0].read(5);
	double err_alpha = parameters_fit_dm_mass[0].read(6);double err_beta =  parameters_fit_dm_mass[0].read(7);

	double mass_[nBin];
	for (int i = 0; i < nBin; i++) mass_[i] = 4.*pi*rho0*pow(r_s,beta)*pow(RLin[0].read(i),(2.-alpha))*pow((RLin[0].read(i)+r_s),(alpha-beta)); 

	array j0, j1, j2, j3;
	j0.init(nBin); j1.init(nBin); j2.init(nBin); j3.init(nBin);

	for (int i = 0; i < nBin; i++) {
		j0.set(i, mass_[i]/rho0);
		j1.set(i, mass_[i] * (1./(RLin[0].read(i)+r_s) * (alpha + beta * (RLin[0].read(i)/r_s))));
		j2.set(i, mass_[i] * log ((RLin[0].read(i)+r_s)/RLin[0].read(i)));
		j3.set(i, mass_[i] * log (r_s/(RLin[0].read(i)+r_s)));
	}
 
	j0 = IntegralTrapezeZeroToPointP(RLin[0], j0);
	j1 = IntegralTrapezeZeroToPointP(RLin[0], j1);
	j2 = IntegralTrapezeZeroToPointP(RLin[0], j2);
	j3 = IntegralTrapezeZeroToPointP(RLin[0], j3);

	for (int i = 0 ; i < nBin ; i++)dm_mass_error[0].set(i, pow( pow(j0.read(i) * err_rho0,2.) + pow (j1.read(i) * err_r_s,2.)+ pow(j2.read(i)*err_alpha,2.)+ pow(j3.read(i)*err_beta,2.),0.5));


}
array err_mass_assume_density_powerLaw(double * fit, array RLin) {
	
	int nBin = RLin.numberBins();
	double rho0 = fit[0];double r_s =  fit[1];
	double alpha = fit[2];double beta =  fit[3];
	
	double err_rho0 = fit[4];double err_r_s = fit[5];
	double err_alpha = fit[6];double err_beta =  fit[7];
	
	double mass_[nBin];
	for (int i = 0; i < nBin; i++) mass_[i] = 4.*pi*rho0*pow(r_s,beta)*pow(RLin.read(i),(2.-alpha))*pow((RLin.read(i)+r_s),(alpha-beta)); 
	
	array j0, j1, j2, j3;
	j0.init(nBin); j1.init(nBin); j2.init(nBin); j3.init(nBin);
	
	for (int i = 0; i < nBin; i++) {
		j0.set(i, mass_[i]/rho0);
		j1.set(i, mass_[i] * (1./(RLin.read(i)+r_s) * (alpha + beta * (RLin.read(i)/r_s))));
		j2.set(i, mass_[i] * log ((RLin.read(i)+r_s)/RLin.read(i)));
		j3.set(i, mass_[i] * log (r_s/(RLin.read(i)+r_s)));
	}
	
	j0 = IntegralTrapezeZeroToPointP(RLin, j0);
	j1 = IntegralTrapezeZeroToPointP(RLin, j1);
	j2 = IntegralTrapezeZeroToPointP(RLin, j2);
	j3 = IntegralTrapezeZeroToPointP(RLin, j3);
	
	array mass_error; mass_error.init(RLin.numberBins());
	for (int i = 0 ; i < nBin ; i++)mass_error.set(i, pow( pow(j0.read(i) * err_rho0,2.) + pow (j1.read(i) * err_r_s,2.)+ pow(j2.read(i)*err_alpha,2.)+ pow(j3.read(i)*err_beta,2.),0.5));
	
	j0.delete_array(); j1.delete_array(); j2.delete_array(); j3.delete_array();
	
	return mass_error;

	
}
array err_density_assume_powerLaw(double * fit, array RLin){
	
	int nBin = RLin.numberBins();
	double rho0 = fit[0];double r_s =  fit[1];
	double alpha = fit[2];double beta =  fit[3];
	
	double err_rho0 = fit[4];double err_r_s = fit[5];
	double err_alpha = fit[6];double err_beta =  fit[7];
	
	double density[nBin];
	for (int i = 0; i < nBin; i++) density[i] = rho0 * pow((r_s + RLin.read(i))/RLin.read(i),alpha) * pow(r_s/(r_s+RLin.read(i)),beta);
	
	array j0, j1, j2, j3;
	j0.init(nBin); j1.init(nBin); j2.init(nBin); j3.init(nBin);
	
	for (int i = 0; i < nBin; i++) {
		
		j0.set(i, density[i] * (1./rho0)); //rho0
		j1.set(i, density[i] * (1./(r_s + RLin.read(i)))* (alpha + ((beta * RLin.read(i))/r_s))); //r_s
		j2.set(i, density[i] * log((RLin.read(i)+r_s)/RLin.read(i))); //alpha
		j3.set(i, density[i] * log(r_s/(RLin.read(i)+r_s))); //beta
	}

	
	array density_error; density_error.init(RLin.numberBins());
	for (int i = 0 ; i < nBin ; i++)density_error.set(i, pow( pow(j0.read(i) * err_rho0,2.) + pow (j1.read(i) * err_r_s,2.)+ pow(j2.read(i)*err_alpha,2.)+ pow(j3.read(i)*err_beta,2.),0.5));
	
	j0.delete_array(); j1.delete_array(); j2.delete_array(); j3.delete_array();
	
	return density_error;
	
}

void cluster::Err_HE_mass_from_data(array Density, array Err_density, array Temperature, array Err_temperature, array Mass_from_HE){
	//change for the values saved at cluster
	mass_HE_error[0] = ArraySQRT(ArrayPlusArray( ArraySQ(ArrayTimesArray(Derivative(Temperature, mass_HE[0]),
															Err_temperature)),
								   ArraySQ(ArrayTimesArray(Derivative(Density, Mass_from_HE),
														   Err_density))));
}
array err_one_variable (array function, array variable, array variable_error){
	array error;
	error.init(function.numberBins());
	error = ArrayTimesArray( ArrayAbs(Derivative(variable,function)),variable_error);
	return error;
}
array err_three_variables(array function, array var1, array var1_error, array var2, array var2_error, array var3, array var3_error){
	
	array error;
	error.init(function.numberBins());
	error = ArraySQRT( ArrayPlusArray(ArraySQ(ArrayTimesArray( ArrayAbs(Derivative(var1,function)),var1_error)),
									  ArrayPlusArray(ArraySQ(ArrayTimesArray( ArrayAbs(Derivative(var2,function)),var2_error)),
													 ArraySQ(ArrayTimesArray( ArrayAbs(Derivative(var3,function)),var3_error)))));
	
	return error;
}

void cluster::Err_dm_beta_from_gamma (){
	dm_beta_error[0] = err_one_variable(dm_beta[0], dm_gamma[0], dm_gamma_error[0]);
}
