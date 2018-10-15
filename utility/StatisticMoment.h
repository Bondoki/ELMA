/*
 * StatisticMoment.h
 *
 *  Created on: 02.12.2014
 *      Author: bondoki
 */

#ifndef STATISTICMOMENT_H_
#define STATISTICMOMENT_H_

#include <cmath>

class StatisticMoment {

private:

	int N; //Anzahl der Elemnte
	double M1; //1.Moment
	double M2; //2.Moment
	double Var; //Varianz
	double sigma; //Standardabweichung
	double Min; // kleinste Element
	double Max; // grooesste Element
	double wert; //Wert weiterreichen

public:
	StatisticMoment() {
		N = 0;
		M1 = 0;
		M2 = 0;
		Var = 0;
		sigma = 0;
		Min = 0;
		Max = 0;
		wert=0;
	}

	/** Brechnung von M1,M2,Var fuer neuen Stichprobenwert
	 *
	 */
	void AddValue(double x) {

		wert = x;
		M1 = M1 * (double) N/(N+1) + x/(N+1);
		M2 = M2 * (double) N/(N+1) + x*x/(N+1);
		N++;
		Var = (M2 - M1*M1)* (double) N/(N-1);
		//Var = (M2 - M1*M1/N)/(N-1);
		sigma = sqrt(Var);


		//if(M1 == Double.NaN)
		//	System.exit(1);

		if ( N == 1)
			Min = x;

		if ( x > Max)
			Max = x;

		if ( x < Min)
			Min = x;
	}

	void clear()
	{
		N = 0;
		M1 = 0;
		M2 = 0;
		Var = 0;
		sigma = 0;
		Min = 0;
		Max = 0;
	}

	const double ReturnM1() const {
		/** Rueckgabe des 1.Momentes
		 */
		//if(M1 == Double.NaN)
		//	System.exit(1);
		return M1;
	}

	const double ReturnM2() const {
		/** Rueckgabe des 2.Momentes
		 */
		return M2;
	}

	double ReturnVar() {
		/** Rueckgabe der Varianz
		 */
		return Var;
	}

	double ReturnSigma() {
		/** Rueckgabe der Standardabeichung
		 */
		return sigma;
	}

	const int ReturnN() const {
		/** Rueckgabe des Stichprobenumfangs
		 */
		return N;
	}

	double ReturnMin() {
		/** Rueckgabe der Standardabeichung
		 */
		return Min;
	}

	double ReturnMax() {
		/** Rueckgabe der Standardabeichung
		 */
		return Max;
	}

	double ReturnWert() {
		/** Rueckgabe des Wertes
		 */
		return wert;
	}
};

#endif /* STATISTICMOMENT_H_ */
