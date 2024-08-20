using System;
using static System.Console;
using static System.Math;

public class root{

	static matrix jacobian(
			Func<vector,vector> f, 
			vector x,
			vector f0 = null
			){ 				//compute the jacobian of a function f
		if(f0==null)f0 = f(x);
		int m = x.size, n = f0.size;
		matrix jac = new matrix(n,m);
		double eps = Pow(2,-26), dx = 1;
		vector x2 = x.copy(), df = new vector(n);
		for(int i=0;i<m;i++){
			dx = eps*(Abs(x2[i])+1e-4);
			x2[i] += dx;
			df = f(x2) - f0;
			jac[i] = df/dx;
			x2[i] = x[i];}
		return jac;
	}//jacobian
	
	public static (vector,vector,int) newton(
			Func<vector,vector> f, 
			vector x0, 
			double eps = 1e-3,
			int max_steps = (int)1e6){
		double lambda = 1; int steps = 0;
		vector x = x0.copy(), Dx = new vector(x0.size), f0 = f(x0), f1 = f0;
                int m = x.size, n = f0.size;
		do{
		steps++;
		matrix J = jacobian(f,x,f0);
		Dx = QRGS.solve(J, -f0);
		lambda = 1;
		f1 = f(x + Dx);
		while(f1.norm() > (1-lambda/2)*f0.norm() && 1 < lambda*1024){
			lambda/=2;
			f1 = f(x+lambda*Dx);
			}
		x+=lambda*Dx;
		f0 = f1;
		}while(f0.norm() >= eps && Dx.norm() >= Pow(2,-26)*x.norm() && steps < max_steps);
		if(steps >= max_steps)Error.WriteLine($"newton: Maximum step of {max_steps} reached");
		return (x,f0,steps);
	}//newton
	
	public static (vector,vector,int) newton_quad_int(
			Func<vector,vector> f, 
			vector x0, 
			double eps = 1e-3,
			int max_steps = (int)1e6){
		double lambda = 1; int steps = 0, _steps = 0;
	        vector x = x0.copy(), Dx = new vector(x0.size), f0 = f(x0), f1 = f0;
	        int m = x.size, n = f0.size;
		double a = 0, b = 0, c = 0;
		do{
                steps++;
		matrix J = jacobian(f,x,f0);
                Dx = QRGS.solve(J, -f0);
		lambda = 1;
	       	f1 = f(x + Dx);	
		_steps = 0;
		c = 0.5*f0.dot(f0);                             //compute quad. interpolation
                b = f0.dot(J*Dx);
		while(f1.norm() > (1-lambda/2)*f0.norm() && _steps<3){		//compute quad. interpolation
				a = (f1.dot(f1)-c)/(lambda*lambda) - b/lambda;
				if(0.1<-b/(2*a) && -b/(2*a)<=1)lambda = -b/(2*a);	//we wish to have lambda in (0,1]
				else lambda/=2;
				f1 = f(x+lambda*Dx);
				_steps++;
				}
		x+=lambda*Dx;
		f0=f1;
                }while(f0.norm() >= eps && Dx.norm() >= Pow(2,-26)*x.norm() && steps <= max_steps);
                if(steps >= max_steps)Error.WriteLine($"newton_quad_int: Maximum step of {max_steps} reached");
                return (x,f0,steps);
	}//newton_quad_int
}//root
