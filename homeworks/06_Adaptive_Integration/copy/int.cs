using System;
using static System.Math;

public static class integrate{
	
	public static (double,double,int) adint(Func<double,double> f, double a, double b,//returns in format (val,err,#eval)
	double delta=0.001, double eps=0.001, double f2=double.NaN, double f3=double.NaN){
		double h=b-a;
		int evals = 0;
		if(Double.IsNaN(f2)){ f2=f(a+2*h/6); f3=f(a+4*h/6); evals+=2;} // first call, no points to reuse
		double f1=f(a+h/6), f4=f(a+5*h/6);
		evals +=2; 				// Add the integrand evaluations
		double Q = (2*f1+f2+f3+2*f4)/6*(b-a); 	// higher order rule
		double q = (  f1+f2+f3+  f4)/4*(b-a); 	// lower order rule
		double err = Abs(Q-q);
		if (err <= delta + eps*Abs(Q))return (Q,err,evals); 	// Accept estimate
			// Bad estimate, do recursion to decrease err.
		(double res1,double err1,int evals1) = adint(f,a,(a+b)/2,delta/Sqrt(2),eps,f1,f2);
		(double res2,double err2,int evals2) = adint(f,(a+b)/2,b,delta/Sqrt(2),eps,f3,f4);   
		evals += evals1 + evals2;
		err = Sqrt(err1*err1 + err2*err2);
		return (res1+res2,err,evals);
	}//adint
	
	public static (double,double,int) transint(Func<double,double> f, double a, double b,
        double delta=0.001, double eps=0.001){ /* applies coord. transform before integrating */
		Func<double,double> ft = delegate(double theta){return f((a+b)/2 + Cos(theta)*(b-a)/2)*Sin(theta)*(b-a)/2;};
		return adint(ft, 0, PI, delta, eps);
	}//transint
	
	public static (double,double,int) integral(Func<double,double> f, double a, double b,//returns in format (val,err,#eval)
        double delta=0.001, double eps=0.001){
		if(double.IsPositiveInfinity(b) && double.IsNegativeInfinity(a)){
			Func<double,double> fs = t=> f(t/(1-t*t))*(1+t*t)/Pow(1-t*t,2);
			return adint(fs,-1,1,delta,eps);}
		if(double.IsPositiveInfinity(b)){
			Func<double,double> fs = t=>f(a + (1-t)/t)/(t*t);
			return adint(fs,0,1,delta,eps);}
		if(double.IsNegativeInfinity(a)){
			Func<double,double> fs = t=>f(b-(1-t)/t)/(t*t);
			return adint(fs,0,1,delta,eps);}
		return adint(f,a,b,delta,eps);
	}//integral

}//integrate
