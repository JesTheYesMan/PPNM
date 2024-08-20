using System;
using static System.Console;
using static System.Math;

public class main{
	public static int Main(){
		Func<vector,double> circ = delegate(vector v){if(v.norm() < 1) return 1; else return 0;};
		Func<vector,double> ball = delegate(vector v){if(v.norm() < 1) return Sqrt(1-v.norm()*v.norm()); else return 0;};
		vector a2d = new vector("-1 -1");
		vector b2d = new vector("1 1");
		vector a3d = new vector("-1 -1 -1");
		vector b3d = new vector("1 1 1");
		int N = 4*(int)1e4;
		double[] exact = new double[3]{PI, 2*PI/3, PI*PI/4};
		double circl = 0, ball2 = 0, ball4 = 0, circ_err = 0, ball2_err = 0, ball4_err = 0;
		for(int i = 500; i < N; i+=500){
			(circl, circ_err) = integrate.quasimc(circ, a2d, b2d, i);
			double circ_eerr = Abs(circl - exact[0]);
                        (ball2, ball2_err) = integrate.quasimc(ball, a2d, b2d, i);
                        double ball2_eerr = Abs(ball2 - exact[1]);
                        (ball4, ball4_err) = integrate.quasimc(ball, a3d, b3d, i);
                        double ball4_eerr = Abs(ball4 - exact[2]);
			Error.WriteLine($"{i} {circ_err} {circ_eerr}  {ball2_err} {ball2_eerr} {ball4_err} {ball4_eerr}");
		}
		Out.WriteLine($"Area of circle, N = {N} exact: {exact[0]}, MC: {circl}, Estimated error: {circ_err}");
                Out.WriteLine($"Volume of half ball, N = {N} exact: {exact[1]}, MC: {ball2}, Estimated error: {ball2_err}");
                Out.WriteLine($"Volume of 4d half ball, N = {N} exact: {exact[2]}, MC: {ball4}, Estimated error: {ball4_err}");
		Out.WriteLine("");
		Func<vector, double> Hf = delegate(vector v){return Pow(PI, -3)/(1-Cos(v[0])*Cos(v[1])*Cos(v[2]));};
		double Hexact = 1.3932039296856768591842462603255;
		vector Ha = new vector("0 0 0");
                vector Hb = new vector(3); for(int i = 0; i < 3; i++) Hb[i] = PI;
		(double Hint, double Herr) = integrate.plainmc(Hf, Ha, Hb, N);
		int M = (int)1e6;
		Out.WriteLine($"Hard integral, N = {M}. MC: {Hint}, est. err. {Herr}, Exact: {Hexact}");
		return 0;
	}//Main
}//main
