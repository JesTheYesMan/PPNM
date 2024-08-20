using System;
using static System.Console;
using static System.Math;

public class main{

	public static int Main(string[] args){
		double rmin = 1e-4, rmax = 8, acc = 0.01, eps = 0.01;
		vector fmin = new vector(2);
		vector fm = new vector(1);
		fmin[0] = rmin*(1-rmin);
		fm[0] = rmin*(1-rmin);
		fmin[1] = 1-2*rmin;
		// M to be minimized without recording wavefunction
		Func<vector,vector> M = delegate(vector v){ 
		Func<double,vector,vector> diff = delegate(double r, vector f){vector w = new vector(2); w[0] = f[1];
                                                                w[1] = (-2.0/r - 2*v[0])*f[0]; return w;};
                vector res = ODE.driver(diff, rmin, fmin, rmax, acc: acc, eps: eps);
		vector tres = new vector(1); tres[0] = res[0]; return tres;};
		// M to be minized, records wave function with eigenvalue
		genlist<double> rs = new genlist<double>();
                genlist<vector> fs = new genlist<vector>();
		
		//Generate and write wavefunction data
		vector Estart = new vector("-1");
		root.newton root = new root.newton(M, Estart);
		vector Es = root.x;
		Func<double,vector,vector> end_diff = delegate(double r, vector f){vector w = new vector(2); w[0] = f[1];
                                                                w[1] = (-2.0/r - 2*Es[0])*f[0]; return w;};
		ODE.driver(end_diff, rmin, fmin, rmax, acc:acc, eps:eps, xlist:rs, ylist:fs);
		for(int i=0;i<rs.size;i++)Error.WriteLine($"{rs[i]} {fs[i][0]}");

		//Start of convergence calculations
		int N = 10;
		vector rmins = new vector(N), rmaxs = new vector(N), accs = new vector(N), epss = new vector(N);
		vector Ermin = new vector(N), Ermax = new vector(N), Eacc = new vector(N), Eeps = new vector(N);
		//construct the values of the parameters
		for(int i=0;i<N;i++){
			rmins[i] = 0.01/Pow(2,i); rmaxs[i] = 6 + i;
			accs[i] = 0.5/Pow(2,i); epss[i] = 0.5/Pow(2,i);}
		//rmin calculations
		for(int i=0;i<N;i++){rmin = rmins[i]; Ermin[i] = new root.newton(M, Estart).x[0];}
		rmin = 1e-4;
		//rmax calculations
		for(int i=0;i<N;i++){rmax = rmaxs[i]; Ermax[i] = new root.newton(M, Estart).x[0];}
                rmax = 8;
		//acc calculations
		for(int i=0;i<N;i++){acc = accs[i]; Eacc[i] = new root.newton(M, Estart).x[0];}
                acc = 0.01;
		//eps calculations
		for(int i=0;i<N;i++){eps = epss[i]; Eeps[i] = new root.newton(M, Estart).x[0];}
                eps = 0.01;

		//Write the data to a file
		for(int i=0;i<N;i++){
			Out.WriteLine($"{rmins[i]} {Ermin[i]} {rmaxs[i]} {Ermax[i]} {accs[i]} {Eacc[i]} {epss[i]} {Eeps[i]}");}
		return 0;
	}//Main
}//main
