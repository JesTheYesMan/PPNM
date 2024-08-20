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

	static double corput(int n, int b){
		double q=0, bk=(double)1.0/b;
		while(n>0){q+=(n%b)*bk; n/=b; bk/=b;}
		return q;
	}//corput
	 
	static void halton(int n, int d, vector x){
		int[] basis = new int[]{2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
		if(d > basis.Length) throw new ArgumentException($"Halton: dimension too large: {d} > {basis.Length}.");
		for(int i = 0; i < d; i++) x[i] = corput(n, basis[i]);
	}//halton
	
	static void latrule(int n, int d, vector x){
		double[] primes = new double[]{2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
		for(int i = 0; i < primes.Length; i++) primes[i] = Sqrt(primes[i]);
		if(d > primes.Length) throw new ArgumentException($"latrule: dimension too large: {d} > {primes.Length}.");
		for(int i = 0; i < d; i++) x[i] = (primes[i] *n) % 1;
	}//latrule

	public static (double,double) plainmc(Func<vector,double> f,vector a,vector b,int N){
        int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
        double sum=0,sum2=0;
	var x=new vector(dim);
	var rnd=new Random();
        for(int i=0;i<N;i++){
                for(int k=0;k<dim;k++)x[k]=a[k]+rnd.NextDouble()*(b[k]-a[k]);
                double fx=f(x); sum+=fx; sum2+=fx*fx;
                }
        double mean=sum/N, sigma=Sqrt(sum2/N-mean*mean);
        var result=(mean*V,sigma*V/Sqrt(N));
        return result;
	}//plainmc
	
	public static (double,double) quasimc(Func<vector,double> f,vector a,vector b,int N){
        int dim=a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
        double sum=0,sum2=0;
        var x=new vector(dim);
	vector halt = new vector(dim);
	vector lat = new vector(dim);
        for(int i=0;i<N;i++){
		halton(i, dim, halt);
		latrule(i, dim, lat);
                for(int k=0;k<dim;k++)x[k]=a[k]+halt[k]*(b[k]-a[k]);
                double fx=f(x); sum+=fx;
		for(int k=0;k<dim;k++)x[k]=a[k]+lat[k]*(b[k]-a[k]);
		fx=f(x); sum2+=fx;
                }
        double mean1 = sum/N, mean2=sum2/N, mean=(mean1+mean2)/2, sigma=Abs(mean1-mean2);
        var result=(mean*V,sigma*V);
        return result;
        }//plainmc
	
	public static (double,double) stratmc(Func<vector,double> f,vector a,vector b,int N, int nmin = 500){
	if(N<nmin) return plainmc(f,a,b,N);
	int dim = a.size; double V=1; for(int i=0;i<dim;i++)V*=b[i]-a[i];
	double sum=0,sum2=0;
        vector x=new vector(dim);
	var rnd = new Random();
	matrix sums = new matrix(2,dim), sums2 = new matrix(2,dim), vars = new matrix(2,dim);
	for(int i=0;i<nmin;i++){						//create nmin random points and add the values in the sums
        for(int k=0;k<dim;k++)x[k]=a[k]+rnd.NextDouble()*(b[k]-a[k]);
		double fx = f(x); sum+=fx; sum2+=fx*fx;
		for(int k=0;k<dim;k++){
		if(x[k] < (b[k] + a[k])/2){sums[0,k]+=fx; sums2[0,k]+=fx*fx;} 	//record the values to the dimension specific memory
		else{sums[1,k]+=fx; sums2[1,k]+=fx*fx;}}
	}
	double mean=sum/nmin, sigma=V*Sqrt(sum2/nmin-mean*mean)/Sqrt(nmin);	//evalute the integral and variance
        int wdim = 0; double maxvar = 0;
	for(int k=0;k<dim;k++){ 						//find the dimension with the largest sub-variance
		vars[0,k] = (sums2[0,k]-sums[0,k]*sums[0,k]/nmin)/nmin;
                vars[1,k] = (sums2[1,k]-sums[1,k]*sums[1,k]/nmin)/nmin;
		if(vars[0,k] > maxvar){maxvar=vars[0,k]; wdim = k;}
                if(vars[1,k] > maxvar){maxvar=vars[1,k]; wdim = k;}}
	vector a_new = a.copy(), b_new = b.copy();				//subdivide the volume
	a_new[wdim] = (a[wdim]+b[wdim])/2; b_new[wdim] = (a[wdim]+b[wdim])/2;

	int N_est = (int)Floor((N-nmin)/(1+vars[0,wdim]/vars[1,wdim]));
	int N_up = Min(N-nmin-2,Max(2,N_est));		//calculate the point distribution, Need at least 2 points in each vol
        int N_down = N-nmin - N_up;
	
	(double int_down, double sigma_down) = stratmc(f, a, b_new, N_down, nmin);
        (double int_up, double sigma_up) = stratmc(f, a_new, b, N_up, nmin);
	double grandint = ((int_down+int_up)*(N-nmin) + V*mean*nmin)/N;	//estimate the grand integral and error
	double grandsigma = Sqrt(Pow(sigma_down*(N-nmin)/N,2) + Pow(sigma_up*(N-nmin)/N,2) + Pow(sigma*nmin/N,2));
	return (grandint,grandsigma);
	}//stratmc

}//integrate
