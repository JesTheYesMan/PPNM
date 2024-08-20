using System;
using static System.Console;
using static System.Math;

public class main{
	
	public static double erf(double x){
	/// single precision error function (Abramowitz and Stegun, from Wikipedia)
	if(x<0) return -erf(-x);
	double[] a={0.254829592,-0.284496736,1.421413741,-1.453152027,1.061405429};
	double t=1/(1+0.3275911*x);
	double sum=t*(a[0]+t*(a[1]+t*(a[2]+t*(a[3]+t*a[4]))));/* the right thing */
	return 1-sum*Exp(-x*x);
	}//erf 

	public static double erfint(double z){ //evaluate erf from integral representation
		if(z < 0) return -erfint(-z);
		if(z > 1){ 
			Func<double,double> f1 = delegate(double t){return Exp(-Pow(z+(1-t)/t,2))/t/t;};
			double res1 = integrate.adint(f1,0,1, 1e-6, 1e-6).Item1;
			return 1-2*res1/Sqrt(PI);}
		Func<double,double> f2 = delegate(double x){return Exp(-x*x);};
		double res2 = integrate.adint(f2,0,z, 1e-7, 1e-7).Item1;
		return 2*res2/Sqrt(PI);
	}//erfint

	public static int Main(){
		//first part of A
		Func<double,double>[] fs = new Func<double,double>[4]{x => Sqrt(x), x => Pow(x,-0.5), 
								x => 4*Sqrt(1-x*x), x => Log(x)/Sqrt(x)};
		double[] Ress = new double[4]{2.0/3.0, 2, PI, -4};
		vector AnalyticalRes = new vector(Ress);
		vector NumRes = new vector(4);
		vector err = new vector(4);
		int[] evals = new int[4];
		WriteLine("Python evaluations are computed using scipy.integrate.quad routine with same requested precision.");
                int[] PyN= {105,231,147,315};
                WriteLine("");
		for(int i = 0; i < 4; i++){
		(NumRes[i],err[i],evals[i]) = integrate.adint(fs[i],0,1,1e-5);
		Out.WriteLine($"{i}: Analytic: {AnalyticalRes[i]}, Numerical: {NumRes[i]}, err: {err[i]}, Evals: {evals[i]}, Python evals: {PyN[i]}");}
		if(vector.approx(NumRes,AnalyticalRes,1e-4)) Out.WriteLine("Numerical integration succesful");
                if(!(vector.approx(NumRes,AnalyticalRes,1e-4))) Out.WriteLine("Numerical integration failed");
		//second part of A
		vector zs = new vector("0 0.02 0.04 0.06 0.08 0.1 0.2 0.3 0.4 0.5 0.6 0.7 0.8 0.9 1 1.1 1.2 1.3 1.4 1.5 1.6 1.7 1.8 1.9 2 2.1 2.2 2.3 2.4 2.5 3.5"); 
		vector erfs = new vector("0 0.022564575 0.045111106 0.067621594 0.090078126 0.112462916 0.222702589 0.328626759 0.428392355 0.520499878 0.603856091 0.677801194 0.742100965 0.796908212 0.842700793 0.880205070 0.910313978 0.934007945 0.952285120 0.966105146 0.976348383 0.983790459 0.989090502 0.992790429 0.995322265 0.997020533 0.998137154 0.998856823 0.999311486 0.999593048 0.999977910 0.999999257");
		int N = 100, M = zs.size;
		for(int i = -N; i < N; i++){
			double z = 3.5*i/(N-1);
			Error.WriteLine($"{z} {erf(z)} {erfint(z)}");
		}
		int NumWin = 0, FunWin = 0;
		double NumErr = 1, FunErr = 1;
		Out.WriteLine("");
		Out.WriteLine($"Comparison of accuracy between Plots approx and numerical integral of Erf:");
		for(int i = 0; i < M; i++){
			NumErr = Abs(erfs[i] - erfint(zs[i])); FunErr = Abs(erfs[i] - erf(zs[i]));
			if(NumErr <= FunErr){ NumWin++;} else FunWin++;
		}
		Out.WriteLine($"Numerical had better approximation {NumWin} times");
                Out.WriteLine($"Plots func had better approximation {FunWin} times");
		Out.WriteLine("");
		for(int i = 0; i < M; i++){
                        NumErr = Abs(erfs[i] - erfint(zs[i])); FunErr = Abs(erfs[i] - erf(zs[i]));
                        Out.WriteLine($"z: {zs[i]} NumErr: {NumErr} FunErr: {FunErr}");
                }
		return 0;
	}//Main
}//main
