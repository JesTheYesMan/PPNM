using System;
using static System.Console;
using static System.Math;

public class main{
	public static int Main(){
		//Three simple functions
		Func<vector,vector> f1 = delegate(vector v){vector w = new vector(1); 
		       					w[0] = v[0]*Exp(-v[0]*v[0]); return w;};
		matrix xs = new matrix("-0.5 0 0.5");
                WriteLine("f = x*e(-x^2). v0 = 0");
                WriteLine("");
                for(int i=0;i<3;i++){
                        vector x = xs[i];
                        root.newton root = new root.newton(f1,x);
			vector xmin = root.x;
			vector fmin = root.f;
			int steps = root.steps;
			int f_eval = root.f_eval;
                        WriteLine($"Numerical results, v0 = {x[0]}:");
                        WriteLine($"v0 = {xmin[0]}. f(v0) = {fmin[0]}");
			WriteLine($"Steps taken: {steps}. Function evals: {f_eval}");
                        WriteLine("");
                }

		Func<vector,vector> f2 = delegate(vector v){vector w = new vector(2); 
							w[0] = 1+v[1]*v[0]; w[1] = 1-v[0]*v[0]; return w;};
                matrix xs2 = new matrix("-0.5 0.001 0.001;0 -0.001 0.5");
                WriteLine("(x,y) -> (1+xy, 1-x^2). v0 = (1,-1) & (-1,1)");
                WriteLine("");
                for(int i=0;i<3;i++){
                        vector x = xs2[i];
                        root.newton root = new root.newton(f2,x);
                        vector xmin = root.x;
                        vector fmin = root.f;
                        int steps = root.steps;
                        int f_eval = root.f_eval;
                        WriteLine($"Numerical results, v0 = ({x[0]}, {x[1]}):");
                        WriteLine($"v0 = ({xmin[0]}, {xmin[1]}). f(v0) = ({fmin[0]}, {fmin[1]})");
                        WriteLine($"Steps taken: {steps}. Function evals: {f_eval}");
                        WriteLine("");
		}

		Func<vector,vector> f3 = delegate(vector v){vector w = new vector(2); 
		       					w[0] = Pow(v[0]-1,2); w[1] = Log(v[0]); return w;};
		matrix xs3 = new matrix("0.5 1.2 2");
                WriteLine("f = ((x-1)^2, ln(x)). x0 = 1");
                WriteLine("");
                for(int i=0;i<3;i++){
                        vector x = xs3[i];
                        root.newton root = new root.newton(f3,x);
                        vector xmin = root.x;
                        vector fmin = root.f;
                        int steps = root.steps;
                        int f_eval = root.f_eval;
                        WriteLine($"Numerical results, x0 = {x[0]}:");
                        WriteLine($"x0 = {xmin[0]}. f(x0) = ({fmin[0]}, {fmin[1]}).");
                        WriteLine($"Steps taken: {steps}. Function evals: {f_eval}");
                        WriteLine("");
                }
		//Rosenbrock function part
		Func<vector,double> Rf = delegate(vector v){return Pow(1-v[0],2)+100*Pow(v[1] - v[0]*v[0],2);};
		Func<vector,vector> dRf = delegate(vector v){vector w = new vector(2);
							w[0] = -2*(1-v[0]) - 400*v[0]*(v[1]-v[0]*v[0]);
							w[1] = 200*(v[1]-v[0]*v[0]); return w;};
		matrix x0s_R = new matrix("2 0 0; 0 0 2");
		WriteLine("f = Rosenbrock function(a=1 b=100). vmin = (1,1), f(vmin) = 0, df(vmin) = 0");
		WriteLine("");
		for(int i=0;i<3;i++){
			vector x0 = x0s_R[i];
                        root.newton root = new root.newton(dRf,x0);
                        vector xmin = root.x;
                        vector dfmin = root.f;
                        int steps = root.steps;
                        int f_eval = root.f_eval;
			double fmin = Rf(xmin);
			WriteLine($"Numerical results, v0 = ({x0[0]}, {x0[1]}):");
			WriteLine($"vmin = ({xmin[0]}, {xmin[1]}). f(vmin) = {fmin}. df(vmin) = ({dfmin[0]}, {dfmin[1]})");
                        WriteLine($"Steps taken: {steps}. Function evals: {f_eval}");
			WriteLine("");
		}

		//Himmelblau function part
                Func<vector,double> Hf = delegate(vector v){;return Pow(v[0]*v[0]+v[1]-11,2)+Pow(v[0] + v[1]*v[1]-7,2);};
                Func<vector,vector> dHf = delegate(vector v){vector w = new vector(2);
                                                        w[0] = 4*v[0]*(v[0]*v[0]+v[1]-11)+2*(v[0]+v[1]*v[1]-7);
                                                        w[1] = 2*(v[0]*v[0]+v[1]-11) + 2*v[1]*(v[0]+v[1]*v[1]-7); 
							return w;};
                matrix x0s_H = new matrix("4 -2 -3 3; 4 2 -2 -2");
                WriteLine("f = Himmelblaus function. vmin = (3,2), (-2.805,3.131), (-3.779,-3.283), (3.584,-1.848)"); 
		WriteLine("f(vmin) = 0, df(vmin) = 0");
                WriteLine("");
                for(int i=0;i<4;i++){
                        vector x0 = x0s_H[i];
                        root.newton root = new root.newton(dHf,x0);
                        vector xmin = root.x;
                        vector dfmin = root.f;
                        int steps = root.steps;
                        int f_eval = root.f_eval;
                        double fmin = Hf(xmin);
                        WriteLine($"Numerical results, v0 = ({x0[0]}, {x0[1]}):");
                        WriteLine($"vmin = ({xmin[0]}, {xmin[1]}). f(vmin) = {fmin}. df(vmin) = ({dfmin[0]}, {dfmin[1]})");
                        WriteLine($"Steps taken: {steps}. Function evals: {f_eval}");
                        WriteLine("");
                        f_eval = 0;
                }
		return 0;
	}//Main
}//main
