/*
 * Tools/Kine
 * namespace Kine for kinematic functions
 * 
 * Author: Tommaso Boschi
 */

#ifndef INTEGRATION_H
#define INTEGRATION_H

#include <limits>

#include "cuba.h"
#include "Tools/asa047.hpp"

namespace Inte
{

	const unsigned int Step = 100;	//high number needed to almost correct integration

	enum MinMax
	{
		True = 1,
		Min = 1,
		Max = -1,
	};

	//integration
	template<class TempClass>
	int Integrand(const int *nDim, const double x[], const int *nComp, double f[], void *UserData)
	{
		TempClass *TempObject = static_cast<TempClass*>(UserData);
		f[0] = TempObject->Function_D(x);
		return 1;
	}

	template<class TempClass>
	double VegasIntegration(TempClass *TempObject, int nDim, int &Trial, int &Fail, double &Error, double &Chi2Prob)
	{
		//input
		double EpsRel = 1.0e-5;		//relative error for each component
		double EpsAbs = 1.0e-6;		//absolute error
		int MinEval = 100;		//minimum number of evaluation
		int MaxEval = 1e5;		//maximum number of evaluation
		int nStart = 20;
		int nIncrease = 10;
		int nBatch = 100;
		void *UserData = TempObject;
		char *state = NULL;
		void *spin = NULL;

		integrand_t IntCast = reinterpret_cast<integrand_t>(&Integrand<TempClass>); 
	
		//output
		double Integral;

		Vegas(nDim, 1, IntCast, UserData, 1, 	//ndim, ncomp, integrand_t, userdata, nvec
		      EpsRel, EpsAbs, 0, 0, 		//epsrel, epsabs, verbosity, seed
		      MinEval, MaxEval, nStart, nIncrease, nBatch,
		      0, state, spin,			//gridno, statefile, spin
		      &Trial, &Fail, &Integral, &Error, &Chi2Prob);
	
		return Integral;
	}

	template<class TempClass>
	double BooleIntegration(TempClass *TempObject)
	{
		double a = 0, b = 0;
		double h = 1.0/Step;
		double Integral = h/90 * 7 * TempObject->Function(0);
		double First = 0;
		for (a = 0; b + 1e-12 < 1.0; a = b)
		{
			b = a + h;

			Integral += First;

			Integral += h/90.0 * 32 * TempObject->Function((3*a +   b) / 4.0);
			Integral += h/90.0 * 12 * TempObject->Function((  a +   b) / 2.0 );
			Integral += h/90.0 * 32 * TempObject->Function((  a + 3*b) / 4.0 );

			First = h/90.0 * 7 * TempObject->Function(b);
			Integral += First;
		}	
	
		//return Integral * h/90.0;
		return Integral;
	}	

	template<class TempClass>
	//double MinLinear(TempClass *TempObject)
	double LinearSolver(TempClass *TempObject, MinMax Sign)
	{
		double MM = TempObject->Function(0);
		double x = 0;

		double h = 1.0/Step;
		for (double a = 0; a < 1.0; a += h)
		{
			double tmp = Sign*TempObject->Function(a);

			if (MM > tmp)
			{
				MM = tmp; 
				x = a;
			}
		}

		return x;
	}	

	template<class TempClass>
	double GoldRatioSolver(TempClass *TempObject, MinMax Sign)
	{
		double GR = (1 + sqrt(5)) / 2.0;

		double S = 0.0;		//start point
		double E = 1.0;		//end point

		while (E - S > 1e-6)
		{
			double A = E - (E - S)/GR;
			double B = S + (E - S)/GR;
			double fA = Sign*TempObject->Function(A);
			double fB = Sign*TempObject->Function(B);

			if (fA < fB)
				E = B;
			else
				S = A;
		}

		return TempObject->Function(Sign*(E - S) <= 0? S : E);
	}

	template<class TempClass>
	double Function(double x[], void *UserData, int Sign)
	{
		TempClass *TempObject = static_cast<TempClass*>(UserData);
		return Sign*TempObject->Function_D(x);
	}

	template<class TempClass>						//num of vars
	double NelMedSolver(TempClass *TempObject, unsigned int n, MinMax Sign)
	//double NelMedSolver(TempClass *TempObject, std::vector<double> &minX, unsigned int n, MinMax Sign)
	{
		//minX.clear();
		void *UserData = TempObject;
		function_t FunCast = reinterpret_cast<function_t>(&Function<TempClass>);

		int icount, ifault, numres;
		int kcount = Step, konvge = 10;
		double reqmin = 1e-5;;
		double *xmin = new double[n];
		double *start = new double[n];	//starting simplex
		double *step = new double[n];	//
		for (unsigned int i = 0; i < n; ++i)
			step[i]  = 0.5;
		start[0] = 0.5;
		start[1] = 0.5;
		start[2] = 0.5;
		start[3] = 0.5;

		double yValue;

		nelmin ( FunCast, n, UserData, Sign, start, xmin, &yValue, 
			 reqmin, step, konvge, kcount, 
			 &icount, &numres, &ifault );

		/*
		if (ifault == 0)
		{
			return Function<TempClass>(xmin, UserData, True);
		}
		else
			return -1.0;
		*/
		return Function<TempClass>(xmin, UserData, True);
	}
}

#endif
