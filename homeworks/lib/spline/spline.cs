using System;
using static System.Console;
using static System.Math;

public class interp{
	public static double linterp(double[] x, double[] y, double z, int i = -1){
		if(i == -1)i = binsearch(x, z); //if statement in case binsearch has been performed previously
		double dx=x[i+1]-x[i]; 
		if(!(dx>0)) throw new Exception("linterp: x-Array improperly ordered");
        	double dy=y[i+1]-y[i];
        	return y[i]+dy/dx*(z-x[i]);
	}//linterp
	 
	public static int binsearch(double[] x, double z){
		if(!(x[0] <= z && x[x.Length - 1] >= z)) throw new ArgumentException($"binsearch: z = {z} outside interval [{x[0]},{x[x.Length-1]}].");
		int i = 0, j = x.Length-1; //endpoints
		while(j-i>1){ 
			int mid=(i+j)/2;
			if(z>x[mid]) i=mid; else j=mid;
		}
		return i;
	}//binsearch
	
	public static double linterpInt(double[] x, double[] y, double z){
		int i = binsearch(x, z);
		double res = 0;
		for(int j = 0; j < i; j++)res += (x[j+1]-x[j])*(y[j+1]+y[j])/2; //add the integrals of each triangle/trapezoid rule
		res += (z-x[i])*(linterp(x,y,z,i) + y[i])/2;
		return res;
	}//linterpInt
	
	public class qspline{
		public vector x,y,b,c;
		public qspline(vector xs,vector ys){
			if(xs.size != ys.size) throw new ArgumentException($"qspline: Data sizes incompatible x: {xs.size}, y:{ys.size}");
			x=xs.copy(); y=ys.copy();
			int n = x.size;
		     	vector p = new vector(n-1), dx = new vector(n-1);
			b = new vector(n-1); 
			c = new vector(n-1);
			for(int i = 0; i < n-1; i++){
				dx[i] = x[i+1] - x[i];
				p[i] = (y[i+1] - y[i])/dx[i];}
			for(int i = 0; i < n-2; i++)c[i+1] =(p[i+1] - p[i] -c[i]*dx[i])/dx[i+1];
			c[n-2] /= 2;
			for(int j = n-3; j >= 0; j--)c[j] = (p[j+1] - p[j] -c[j+1]*dx[j+1])/dx[j];
			for(int i = 0; i < n-1; i++)b[i] = p[i] - c[i]*dx[i];
		}//constructor
		
		public double evaluate(double z){
			int i = binsearch(x, z);
			double res = y[i] + b[i]*(z - x[i]) + c[i]*(z-x[i])*(z - x[i]);
			return res;
		}//evaluate

		public double derivative(double z){
			int i = binsearch(x, z);
			double res = b[i] + 2*c[i]*(z-x[i]); 
			return res;
		}//derivative

		public double integral(double z){
			int i=binsearch(x,z);
			double res = 0;
			for(int j = 0; j < i; j++){
				res += (x[j+1] - x[j])*y[j] + b[j]*Pow(x[j+1]-x[j],2)/2 + c[i]*Pow(x[j+1] - x[j],3)/3;}
			res+= (z - x[i])*y[i] + b[i]*Pow(z-x[i],2)/2 + c[i]*Pow(z - x[i],3)/3;
			return res;	
		}//integral
	}//qspline
	 
	public class cspline{
		public vector x,y,b,c,d;
		public cspline(vector xs, vector ys){
			if(xs.size != ys.size) throw new ArgumentException($"cspline: Data sizes incompatible x: {xs.size}, y:{ys.size}");
			x = xs.copy(); y = ys.copy();
			int n = x.size;
			//construct all relevant vectors and matrices
			vector h = new vector(n-1), p = new vector(n-1), B = new vector(n);
			b = new vector(n); c = new vector(n-1); d = new vector(n-1);
			matrix A = new matrix(n);
			for(int i = 0; i < n-1; i++){
				h[i] = x[i+1]-x[i];
				p[i] = (y[i+1] - y[i])/h[i];}
			A[0,0] = 2; A[n-1,n-1] = 2; A[0,1] = 1; B[0] = 3*p[0]; B[n-1] = 3*p[n-2]; A[n-1,n-2] = 1;
			for(int i = 0; i < n-2; i++){
				A[i+1,i+1] = 2*h[i]/h[i+1] + 2;
				A[i+1,i+2] = h[i]/h[i+1];
				A[i+1,i] = 1;
				B[i+1] = 3*(p[i]+p[i+1]*h[i]/h[i+1]);}
			//Now solve for b using tridiag-solver
			b = linsol.TriDiagSol(A,B);
			//compute c and d using eq. 19-20 in the book
			c[0] = 0;
			for(int i = 0; i < n-2; i++){
				c[i+1] =(-2*b[i+1] -b[i+2]+3*p[i+1])/h[i+1];
				d[i] = (b[i] +b[i+1] - 2*p[i])/(h[i]*h[i]);}
			d[n-2] = -c[n-2]/(3*h[n-2]);
		}//constructor
		
		public double evaluate(double z){
                        int i = binsearch(x, z);
                        double res = y[i] + b[i]*(z-x[i]) + c[i]*Pow(z-x[i],2) + d[i]*Pow(z-x[i],3);
                        return res;
                }//evaluate

		public double derivative(double z){
                        int i = binsearch(x, z);
                        double res = b[i] + 2*c[i]*(z-x[i]) + 3*d[i]*(z-x[i])*(z-x[i]);
                        return res;
                }//derivative

                public double integral(double z){
                        int i=binsearch(x,z);
                        double res = 0;
			for(int j = 0; j < i; j++){ //split terms for ease of reading
				double T12 = (x[j+1] - x[j])*y[j] + b[j]*Pow(x[j+1]-x[j],2)/2;
				double T34 = c[j]*Pow(x[j+1] - x[j],3)/3 + d[j]*Pow(x[j+1]-x[j],4)/4;
                                res += T12 + T34;}
                        res+= (z - x[i])*y[i] + b[i]*Pow(z-x[i],2)/2 + c[i]*Pow(z - x[i],3)/3 + d[i]*Pow(z-x[i],4)/4;
                        return res;
                }//integral
	}

	public class vcspline{ // spline vectors vecy(x). Simplest implementation, not necessarily best
		public vector x;
		public vector[] y;
		public cspline[] splines;
		public int n, dim;
		public vcspline(vector xs, vector[] ys){
			if(ys[0].size!=xs.size) throw new ArgumentException($"vcspline: x({xs.size}) and y({ys[0].size}) must have same size.");
			n = xs.size; dim = ys.Length;
			x = xs.copy();
			y = new vector[dim];
			splines = new cspline[dim];
			for(int i=0;i<dim;i++){
				y[i] = ys[i].copy();
				splines[i] = new cspline(x,y[i]);
			}
		}//constructor

		public vector evaluate(double z){
			vector res = new vector(n);
			for(int i=0;i<dim;i++)res[i]=splines[i].evaluate(z);
			return res;
		}//evaluate
		
		public vector derivative(double z){
                        vector res = new vector(n);
                        for(int i=0;i<dim;i++)res[i]=splines[i].derivative(z);
                        return res;
                }//derivative

		public vector integral(double z){
                        vector res = new vector(n);
                        for(int i=0;i<dim;i++)res[i]=splines[i].evaluate(z);
                        return res;
                }//integral
	}//vcspline 
}//spline
