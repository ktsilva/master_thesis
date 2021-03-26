/*
 *  arrayMathematics.h
 *  
 *
 *  Created by Catarina Silva Fernandes 2010.
 *  No copyrights! Go for it and use it freely!
 *
 */


/* "Mathematics is the language with which God has written the universe." Galileu
 * Creates all necessary operations needed to deal with arrays
 * Nothing fancy: copy an array to array, sums arrays, array with double,
 * makes all sort of operations, sine, cosine...
 * For fancy: you can find derivative and integral
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
#include "gsl/gsl_rng.h"
#include "rng.h"

using namespace std;
string choice; //let us put this somewhere else

double rdmNumber (double min, double max) {
	
	return cat_rng_uniform_double_min_max(min, max);
}

/*int intRdmNumber(int min, int max){
	
	return cat_rng_uniform_int_min_max( min,  max);
}*/

double Likelihood(array A, array B){
	
	int a = A.numberBins();
	
	double chiSq = 0.;
	
	for (int i = 0; i < a; i++){
		chiSq+=pow( (A.read(i) - B.read(i))/(A.read(i) + B.read(i)),2);
		printf("ChiSq = %g\n", chiSq);
	}
	
	
	return chiSq;
	
	
}


//Sets elements to an array
array CopyArray(array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	for (int i = 0; i < a; i++) {
		C.set(i, A.read(i));
	}
	return C;
}
array CopyArrayMagnifiedDouble (array A, double magnify) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	for (int i = 0; i < a; i++) {
		C.set(i, magnify * A.read(i));
	}
	return C;
}
array CopyArrayMagnifiedInt (array A, int magnify) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	for (int i = 0; i < a; i++) {
		C.set(i, 1.* magnify * A.read(i));
	}
	return C;
}
array SetAllToDouble (int numberBins, double b) {
	int a = numberBins;
	
	array C;
	C.init(a);
	
	for (int i = 0; i < a; i++) {
		C.set(i, b);
	}
	return C;
}
array SetAllToInt (int numberBins, int b) {
	int a = numberBins;
	
	array C;
	C.init(a);
	C = SetAllToDouble(a, 1.*b);
	
	return C;
}
array SetAllToZero (int numberBins) {
	int a = numberBins;
	
	array C;
	C.init(a);
	
	C = SetAllToDouble(a, 0.);
	
	return C;
}
array EraseArray(int numberBins) {
	int a = numberBins;
	
	array C;
	C.init(a);
	
	C = SetAllToDouble(a, 0.);
	
	return C;
}
array ArrayAbs(array A){
	
	array C;
	C.init(A.numberBins());
	
	for (int i = 0; i < A.numberBins(); i++) C.set(i, abs(A.read(i)));
	
	return C;
}
//write to File
string datetime_to_string(const tm& time, const char* format) 
{
    stringstream datetime;
	
    // retrieve the time_put facet installed in the stream
    const time_put<char>& writer = 
	use_facet< time_put<char> >(datetime.getloc());
	
    int len = strlen(format);
	
    //formats the contents of the tm time into the output stream datetime
    if (writer.put(datetime, datetime, ' ', 
				   &time, format, format + len).failed( )) 
	{
        printf("formatting date time failed!");
	}
	
    return datetime.str();
}
void writeArrayToFile (array xAxis, array yAxis, string nameFile) {
	
	int x = xAxis.numberBins();
	
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
	
	for (int i = 0; i < x; i++) {
		
		outFile << xAxis.read(i) << " " << yAxis.read(i) << endl;
		
	}
	
	outFile.close();
}

void writeDoubleToFile (double *xAxis, double *yAxis, string nameFile, int Nslots) {
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
	
	for (int i = 0; i < Nslots; i++) {
		
		outFile << xAxis[i] << " " << yAxis[i] << endl;
		
	}
	
	outFile.close();
}

void write_array_and_error_to_file(array xAxis, array yAxis, array yError, string nameFile) {
	
	int x = xAxis.numberBins();
	
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
	
	for (int i = 0; i < x; i++) {
		
		outFile << xAxis.read(i) << " " << yAxis.read(i) << " " << yError.read(i) << endl;
		
	}
	
	outFile.close();
}

//Element+Element
array ArrayPlusArray (array A, array B) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = A.read(i) + B.read(i);
		C.set(i, value);
	}
	
	return C;
}
array ArrayPlusDouble (array A, double b) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = A.read(i) + b;
		C.set(i, value);
	}
	return C;
}
array ArrayPlusInt (array A, int b) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	C = ArrayPlusDouble(A, 1.*b);
	
	return C;
	
	
}

//Element-Element
array ArrayMinusArray (array A, array B) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = A.read(i) - B.read(i);
		C.set(i, value);
	}
	
	return C;
}
array ArrayMinusDouble (array A, double b) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = A.read(i) - b;
		C.set(i, value);
	}
	
	return C;
}
array ArrayMinusInt (array A, int b) {
	
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	C = ArrayMinusDouble (A, 1.*b);
	
	return C;
}
array DoubleMinusArray (double b, array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = b - A.read(i);
		C.set(i, value);
	}
	
	return C;
}
array IntMinusArray (int b, array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	C = DoubleMinusArray (1.*b, A);
	
	return C;
}

// Element*Element
array ArrayTimesArray (array A, array B) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = A.read(i) * B.read(i);
		C.set(i, value);
	}
	
	return C;
}
array DoubleTimesArray (double b, array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = A.read(i) * b;
		C.set(i, value);
	}
	
	return C;
	return C;	
}
array IntTimesArray (int b, array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	C = DoubleTimesArray (1.*b, A);
	
	return C;
}

// Nominator/Denominator
array ArrayDividesArray (array nominator, array denominator) {
	int a = nominator.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = nominator.read(i) / denominator.read(i);
		C.set(i, value);
	}
	
	return C;
}
array ArrayDividesDouble (array nominator, double denominator) {
	int a = nominator.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = nominator.read(i) / denominator;
		C.set(i, value);
	}
	
	return C;
}
array ArrayDividesInt (array nominator, int denominator) {
	int a = nominator.numberBins();
	
	array C;
	C.init(a);
	
	C = ArrayDividesDouble (nominator, 1. * denominator);
	
	return C;
}
array DoubleDividesArray (double nominator, array denominator) {
	int a = denominator.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = nominator / denominator.read(i);
		C.set(i, value);
	}
	
	return C;
}
array IntDividesArray (int nominator, array denominator) {
	int a = denominator.numberBins();
	
	array C;
	C.init(a);
	
	C = DoubleDividesArray ( 1. * nominator, denominator);
	
	return C;
}

//Mathematical functions for array A
array ArraySQ (array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = A.read(i) * A.read(i);
		C.set(i, value);
	}
	
	return C;
}
array ArraySQRT (array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = pow(A.read(i),0.5);
		C.set(i, value);
	}
	
	return C;
}
array ArrayPowerDouble (array A, double power) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = pow(A.read(i), power);
		C.set(i, value);
	}
	
	return C;
}
array ArrayPowerInt (array A, int power) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	C = ArrayPowerDouble (A, 1. * power);
	
	return C;
}
array doublePowerArray (double b, array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = pow(b,A.read(i));
		C.set(i, value);
	}
	
	return C;
}
array tenPowerArray (array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	
	C = doublePowerArray (10., A);
	
	return C;
}
array expPowerArray (array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = exp(A.read(i));
		C.set(i, value);
	}
	
	return C;
}
array Log10Array (array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = log10(A.read(i));
		C.set(i, value);
	}
	
	return C;
}
array LnArray (array A) {
	int a = A.numberBins();
	
	array C;
	C.init(a);
	double value = 0.;
	
	for (int i = 0; i < a; i++) {
		value = log(A.read(i));
		C.set(i, value);
	}
	
	return C;
}
array SineArray (array A) {
	int a = A.numberBins();
	
	array C;	
	C.init(a);
	
	double value = 0.; 
	
	for (int i = 0; i < a;i++) {
		value = sin(A.read(i));
		C.set(i, value);
	}
	
	return C;
}
array CosineArray (array A) {
	int a = A.numberBins();
	
	array C;	
	C.init(a);
	
	double value = 0.; 
	
	for (int i = 0; i < a;i++) {
		value = cos(A.read(i));
		C.set(i, value);
	}
	return C;
}

//Retuns an array in specific scale given min and max value
array setLog10scale (double scaleMin, double scaleMax, int numberBins) {
	int a = numberBins;
	
	array C;
	C.init(a);
	
	double delta = (log10(scaleMax) - log10(scaleMin))/a;
	double exponent = log10(scaleMin);
	
	for (int i = 0; i < numberBins; i++) {
		exponent += delta;	
		C.set (i, exponent);
	}
	
	return C;
	
}
array setLnScale (double scaleMin, double scaleMax, int numberBins) {
	int a = numberBins;
	
	array C;
	C.init(a);
	
	double delta = (log(scaleMax) - log(scaleMin))/a;
	double exponent = log(scaleMin);
	
	for (int i = 0; i < numberBins; i++) {
		exponent += delta;	
		C.set (i, exponent);
	}
	
	return C;
	
}
array SetLinScale (int numberBins, double maxValue) {	// set(0, delta) 
	int a = numberBins;
	
	array C;
	C.init(numberBins);
	double delta = maxValue/numberBins;
	
	for (int i = 0; i < a; i++) {
		C.set(i, (i+1) * delta);
	}
	
	return C;
}

//Take a closer look into the array
void coutOneArrayToScreen (array Out) {
	int o = Out.numberBins();
	
	for (int i = 0; i < o; i++) {
		cout << i << " " << Out.read(i) << endl;
	}
}
void coutTwoArrayToScreen (array Out1, array Out2) {
	int o = Out1.numberBins();
	
	for (int i = 0; i < o; i++) {
		cout << i << " " << Out1.read(i) << " " << Out2.read(i) << endl;
	}
}
void coutThreeArrayToScreen (array Out1, array Out2, array Out3) {
	int o = Out1.numberBins();
	
	for (int i = 0; i < o; i++) {
		cout << i << " " << Out1.read(i) << " " << Out2.read(i) << " " << Out3.read(i) << endl;
	}
}

//Sets ground for derivative and integral calculus
double SlopePointPandQ (array xAxis, array yAxis, int PointP, int PointQ) {
	int i = PointP;
	int j = PointQ;
	
	double slope;
	slope = (yAxis.read(j)-yAxis.read(i)) / (xAxis.read(j) - xAxis.read(i));
	
	return slope;
	
}
double DerivativePointP (array xAxis, array yAxis, int PointP) {
	int x = xAxis.numberBins();
	int i = PointP;
	
	double derivative, yy, xx, yLeft, xLeft, yRight, xRight;
	
	if (i != 0 && i != x-1) {
		
		yRight = yAxis.read(i+1);
		yLeft = yAxis.read(i-1);
		yy = yAxis.read(i);
		
		xRight = xAxis.read(i+1);
		xLeft = xAxis.read(i-1);
		xx = xAxis.read(i);
		
		derivative = 0.5 * (((yy - yLeft)/(xx - xLeft))+((yRight - yy)/(xRight - xx)));
	}
	
	if (i == 0) { //attention was given
		yRight = yAxis.read(i+1);
		yy = yAxis.read(i);
		
		xRight = xAxis.read(i+1);
		xx = xAxis.read(i);
		
		derivative = (yRight - yy)/(xRight - xx);
	}
	
	if (i == x-1) {	//attention was given
		yy = yAxis.read(i);
		yLeft = yAxis.read(i-1);
		
		xx = xAxis.read(i);
		xLeft = xAxis.read(i-1);
		
		derivative = (yy - yLeft)/(xx - xLeft);
	}
	
	return derivative;
}
array Derivative (array xAxis, array yAxis) {
	int x = xAxis.numberBins();
	
	array C;
	C.init(x);
	
	for (int i = 0; i < x; i++)  {
		
		C.set(i, DerivativePointP (xAxis, yAxis, i));
	}
	
	return C;
}

array Gamma(array xAxis, array yAxis){

	int x = xAxis.numberBins();
	
	array C;
	C.init(x);
	
	C = Derivative(Log10Array(xAxis), Log10Array(yAxis));
	
	return C;
	
}

array DerivativeAndCorrect (array xAxis, array yAxis) {
	int x = xAxis.numberBins();
	
	array C;
	C.init(x);
	
	for (int i = 0; i < x; i++)  {
		
		C.set(i, DerivativePointP (xAxis, yAxis, i));
	}
	
	/*	CORRECTION
	 */
	
	double slopeLeft, slopeRight, slope, valueTrue, slopeRdm;
	cout << "correct first values [y/n]" << endl;
	cin >> choice;
	double error;
	cin >> error;
	
	if (choice == "y" ) {
		for (int i = x*0.5; i>-1; i--) {
			
			slopeRight = SlopePointPandQ(xAxis, C, i+1, i+3);
			slope = SlopePointPandQ (xAxis, C, i,i+2);
			
			
			
			if (slopeRight > slope && (slopeRight-slope)/slopeRight > error ) {
				slopeRdm = slopeRight * rdmNumber (1.-error, 1.);
				valueTrue = C.read(i+2) - slopeRdm * (xAxis.read(i+2) - xAxis.read(i));
				
				C.set (i, valueTrue);
				
			}
			
			if (slopeRight < slope && (slopeRight-slope)/slopeRight < -1.*error){
				slopeRdm = slopeRight * rdmNumber (1., 1.+error);
				valueTrue = C.read(i+2) - slopeRdm * (xAxis.read(i+2) - xAxis.read(i));
				
				C.set (i, valueTrue);
				
				
			}
			
		}
	}
	for (int i = x*0.5; i<x; i++) {
		
		/*	slopeLeft = Derivative (xAxis, C, i+1);
		 slope = Derivative (xAxis, C, i);
		 slopeRight = Derivative (xAxis, C, i-1);
		 */	
		slopeLeft = SlopePointPandQ(xAxis, C, i-3, i-1);
		slope = SlopePointPandQ (xAxis, C, i-2,i);
		
		if (slopeLeft > slope && (slopeLeft-slope)/slopeLeft > error) {
			slopeRdm = slopeLeft * rdmNumber (1.-error, 1.);
			valueTrue = slopeRdm * (xAxis.read(i) - xAxis.read(i-2)) + C.read(i-2);
			C.set (i, valueTrue); 
			
		}
		if (slopeLeft < slope && (slopeLeft-slope)/slopeLeft < -1.*error){
			slopeRdm = slopeLeft * rdmNumber (1., 1.+error);
			valueTrue = slopeRdm * (xAxis.read(i) - xAxis.read(i-2)) + C.read(i-2);
			C.set (i, valueTrue); 
		}
	}
	
	
	
	return C;
}
array IntegralTrapezeZeroToPointP(array xAxis, array yAxis) { //all the points
	int x = xAxis.numberBins();
	
	array C;
	C.init(x);
	
	double integral = 0.;
	double yInit, yFinal, xInit, xFinal;
	
	for (int i = 0; i < x; i++)  {
		if (i == 0) { //attention was given
			yFinal = yAxis.read(i);
			yInit = yFinal;
			
			xFinal = xAxis.read(i);
			xInit = 0.;
		}
		
		else {
			yFinal = yAxis.read(i);
			yInit = yAxis.read(i-1);
			
			xFinal = xAxis.read(i);
			xInit = xAxis.read(i-1);
		} 
		
		integral += 0.5 * (xFinal - xInit) * (yFinal + yInit);
		C.set(i, integral);
	}
	return C;
}
double * integralTrapezeZeroToPointP(int nBins, double * xAxis, double * yAxis){
	int x = nBins;
	
	double * c;
	c = new double [x];
	
	double integral = 0.;
	double yInit, yFinal, xInit, xFinal;
	
	for (int i = 0; i < x; i++)  {
		if (i == 0) { //attention was given
			yFinal = yAxis[i];
			yInit = yFinal;
			
			xFinal = xAxis[i];
			xInit = 0.;
		}
		
		else {
			yFinal = yAxis[i];
			yInit = yAxis[i-1];
			
			xFinal = xAxis[i];
			xInit = xAxis[i-1];
		} 
		
		integral += 0.5 * (xFinal - xInit) * (yFinal + yInit);
		c[i] = integral;
	}
	return c;
	
}
array IntegralTrapezePointPtoInfty (array xAxis, array yAxis) {
	int x = xAxis.numberBins();
	
	array C;
	C.init(x);
	
	double integral = 0.;
	double xF, yF;
	double xI, yI;
	
	for (int i = x - 1; i > 0; i--) {
		xF = xAxis.read(i);
		yF = yAxis.read(i);
		xI = xAxis.read(i-1);
		yI = yAxis.read(i-1) ;
		
		integral += 0.5 * (xF - xI) * (yF + yI);
		
		C.set(i, integral);
	}
	C.set(0, C.read(1));
	
	
	return C;	
}

//Returns y' from the followin general Eq.: dlogy' = dlogy + f_x dlogx
array YyFromdLogYydLogX (array X, array F_x, double Yy0, array Y) {
	
	int x = X.numberBins();
	
	array Yy, integral;
	Yy.init(x);
	integral.init(x);
	
	integral = IntegralTrapezeZeroToPointP (LnArray(X), F_x);
	
	Yy.set(0, Yy0);
	
	double value = 0.;
	for (int i = 1; i < x; i++) {
		value = (Yy0/Y.read(0)) * Y.read(i)* exp(integral.read(i));
		Yy.set(i, value);
	}
	
	integral.delete_array();
	
	return Yy;
	
}

//Can be used for numerical purposes if the extremes look too crazy - needs to be improved
array correct (array xAxis, array yAxis, double aScale){
	
	int x = xAxis.numberBins();
	
	array C;
	C.init(x);
	
	double derivativeR, derivativeL;
	double value =0.;
	
	for (int i = x-1; i > -1; i--) {
		if ( xAxis.read(i) < aScale) {
			derivativeR = DerivativePointP(xAxis, yAxis, i+1);
			derivativeL = DerivativePointP(xAxis, yAxis, i);
			if (abs((derivativeL-derivativeR)/derivativeR) > 0.01) {
				value = C.read(i+1) - (derivativeR * (xAxis.read(i+1) - xAxis.read(i)));
				C.set(i, value);
			}
			else {
				C.set(i, yAxis.read(i));
			}
		}
		else{
			C.set(i, yAxis.read(i));
		}
	}
	
	
	return C;
}
array correctInitFinal(array yAxis, int RMIN, int RMAX) {
	
	int y = yAxis.numberBins();
	double valueMin = yAxis.read(RMIN);
	double valueMax = yAxis.read(RMAX);
	for (int i = 0; i < RMIN; i++) {
		yAxis.set(i,valueMin);
	}
	for (int i = RMAX+1; i < y; i++) {
		yAxis.set(i,valueMax);
	}
	
	return yAxis;
	
}