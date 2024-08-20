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
	}
	public static vector driver(
	    Func<double, vector, vector> f,
	    double a,
	    vector ya,
	    double b,
	    double h = 0.01,
	    double hmax = Double.NaN,
	    double acc = 0.01,
	    double eps = 0.01,
	    int nmax = 9999,
	    genlist<double> xlist = null,
	    genlist<vector> ylist = null
	)
	{
	    if (a > b) throw new ArgumentException("driver: a>b");
	
	    double x = a;
	    vector y = ya.copy();
	    vector tol = new vector(y.size);
	
	    if (xlist != null) xlist.add(x);
	    if (ylist != null) ylist.add(y);
	
	    int steps = 0;
	    do
	    {
	        if (x >= b) return y;
	        if (x + h > b) h = b - x;
	
	        (var yh, var err) = rkstep12(f, x, y, h);
	
	        for (int i = 0; i < y.size; i++)
	            tol[i] = (acc + eps * Abs(yh[i])) * Sqrt(h / (b - a));
	
	        bool ok = true;
	        for (int i = 0; i < y.size; i++)
	            if (!(err[i] < tol[i])) ok = false;
	
	        if (ok)
	        {
	            x += h;
	            y = yh;
	
	            if (xlist != null) xlist.add(x);
	            if (ylist != null) ylist.add(y);
	        }
	
	        double factor = tol[0] / Abs(err[0]);
	        for (int i = 1; i < y.size; i++)
	            factor = Min(factor, tol[i] / Abs(err[i]));
	
	        h *= Min(Pow(factor, 0.25) * 0.95, 2);
	        if (!Double.IsNaN(hmax) && h > hmax) h = hmax;
	
	        steps++;
	    } while (steps <= nmax);
	
	    if (xlist != null) return y;
	    throw new ArgumentException($"Driver: Did not finish within allotted steps, time reached {x}");
	}
	
	public static interp.vcspline driver_interp(
	    Func<double, vector, vector> f,
	    double a,
	    vector ya,
	    double b,
	    double h = 0.01,
	    double hmax = Double.NaN,
	    double acc = 0.01,
	    double eps = 0.01,
	    int nmax = 9999
	)
	{
	    var x = new genlist<double>();
	    var y = new genlist<vector>();
	    int dim = ya.size;	
	    driver(f, a, ya, b, h, hmax, acc, eps, nmax, x, y);	
	    var x_data = new vector(x.size);
	    var y_data = new vector[dim];
	    for (int i = 0; i < x.size; i++)
	        x_data[i] = x[i];
	
	    for (int j = 0; j < dim; j++)
	    {
	        y_data[j] = new vector(x.size);
	        for (int i = 0; i < x.size; i++)
	            y_data[j][i] = y[i][j];
	    }	
	    return new interp.vcspline(x_data, y_data);
	}
}
