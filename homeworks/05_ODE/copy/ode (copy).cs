using System;
using static System.Math;

public static class ODE{
	public static (vector, vector) rkstep12(
			Func<double,vector,vector> f, 	/* dy/dx = f(x,y) <- */
			double x, 			/* current value of x */
			vector y, 			/* current value of y */
			double h			/* step size */
			){
			vector k0 = f(x,y);              /* embedded lower order formula (Euler) */
		vector k1 = f(x+h/2,y+k0*(h/2)); /* higher order formula (midpoint) */
		vector yh = y+k1*h;              /* y(x+h) estimate */
		vector er = (k1-k0)*h;           /* error estimate */
		return (yh,er);
	}//rkstep12
	
	public static vector driver(
		Func<double,vector,vector> f,/* the f from dy/dx=f(x,y) */
		double a,                    /* the start-point a */
		vector ya,                   /* y(a) */
		double b,                    /* the end-point of the integration */
		double h=0.01,               /* initial step-size */
		double hmax = Double.NaN,    /* maximal allowed stepsize */
		double acc=0.01,             /* absolute accuracy goal */
		double eps=0.01,             /* relative accuracy goal */
		int nmax = 9999,	     /* Max allowed steps */
		genlist<double> xlist=null,  /* Initialized x list if path needs to be recorded*/
                genlist<vector> ylist=null   /* Initialized y list if path needs to be recorded*/
		){
		if(a>b) throw new ArgumentException("driver: a>b");
		double x=a; vector y=ya.copy(), tol = new vector(y.size);
		if(xlist!=null)xlist.add(x);
		if(ylist!=null)ylist.add(y);
		int steps = 0;
		do{
	        if(x>=b) return y; /* job done */
        	if(x+h>b) h=b-x;   /* last step should end at b */
        	(var yh,var err) = rkstep12(f,x,y,h);
        	for(int i=0;i<y.size;i++)tol[i]=(acc+eps*Abs(yh[i]))*Sqrt(h/(b-a)); /* Evaluate the tolerances*/
                bool ok=true;
                for(int i=0;i<y.size;i++)if(!(err[i]<tol[i])) ok=false; /* check whether to accept step */
                if(ok){ 						/* step accepted */
			x+=h; y=yh;
		   	if(xlist!=null)xlist.add(x);			/* record path if desired */
			if(ylist!=null)ylist.add(y);} 
                double factor = tol[0]/Abs(err[0]);
                for(int i=1;i<y.size;i++) factor=Min(factor,tol[i]/Abs(err[i])); /* figure out new step size*/
                h *= Min( Pow(factor,0.25)*0.95 ,2);
		if(hmax!=Double.NaN && h > hmax)h=hmax;
		steps++;
        	}while(steps<=nmax);
		if(xlist!=null) return y;
		throw new ArgumentException($"Driver: Did not finish within alloted steps, time reached {x}");
	}//driver
	
	public static interp.vcspline driver_interp( /*returns a vector cubic spline of the solution, room for optimization */
                Func<double,vector,vector> f,/* the f from dy/dx=f(x,y) */
                double a,                    /* the start-point a */
                vector ya,                   /* y(a) */
                double b,                    /* the end-point of the integration */
                double h=0.01,               /* initial step-size */
                double hmax = Double.NaN,    /* maximal allowed stepsize */
                double acc=0.01,             /* absolute accuracy goal */
                double eps=0.01,             /* relative accuracy goal */
                int nmax = 9999              /* Max allowed steps */
                ){
		genlist<double> x = new genlist<double>();
		genlist<vector> y = new genlist<vector>();
	        int dim = ya.size;
		driver(f,a,ya,b,h,hmax,acc,eps,nmax,x,y);
		vector x_data = new vector(x.size);
		vector[] y_data = new vector[dim];
		for(int i=0;i<x.size;i++)x_data[i] = x[i];
			for(int j=0;j<dim;j++){
                                y_data[j] = new vector(x.size);
				for(int i=0;i<x.size;i++)y_data[j][i] = y[i][j];
			}
		interp.vcspline res = new interp.vcspline(x_data, y_data);
		return res;
	}//driver_interp
}//ODE
