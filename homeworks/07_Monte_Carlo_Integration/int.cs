using System;
using static System.Math;

public class integrate{

    static double corput(int n, int b){
        double q = 0, bk = 1.0 / b;
        while(n > 0){
            q += (n % b) * bk; 
            n /= b; 
            bk /= b;
        }
        return q;
    }
    
    static void halton(int n, int d, vector x){
        int[] basis = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
        if(d > basis.Length) throw new ArgumentException($"Halton: dimension too large: {d} > {basis.Length}.");
        for(int i = 0; i < d; i++) x[i] = corput(n, basis[i]);
    }
    
    static void latrule(int n, int d, vector x){
        double[] primes = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};
        for(int i = 0; i < primes.Length; i++) primes[i] = Sqrt(primes[i]);
        if(d > primes.Length) throw new ArgumentException($"latrule: dimension too large: {d} > {primes.Length}.");
        for(int i = 0; i < d; i++) x[i] = (primes[i] * n) % 1;
    }

    public static (double, double) plainmc(Func<vector, double> f, vector a, vector b, int N){
        int dim = a.size; 
        double V = 1; 
        for(int i = 0; i < dim; i++) V *= b[i] - a[i];
        double sum = 0, sum2 = 0;
        var x = new vector(dim);
        var rnd = new Random();
        for(int i = 0; i < N; i++){
            for(int k = 0; k < dim; k++) x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
            double fx = f(x); 
            sum += fx; 
            sum2 += fx * fx;
        }
        double mean = sum / N, sigma = Sqrt(sum2 / N - mean * mean);
        var result = (mean * V, sigma * V / Sqrt(N));
        return result;
    }
    
    public static (double, double) quasimc(Func<vector, double> f, vector a, vector b, int N){
        int dim = a.size; 
        double V = 1; 
        for(int i = 0; i < dim; i++) V *= b[i] - a[i];
        double sum1 = 0, sum2 = 0;
        vector x = new vector(dim), halt = new vector(dim), lat = new vector(dim);
        for(int i = 0; i < N; i++){
            halton(i, dim, halt);
            latrule(i, dim, lat);
            for(int k = 0; k < dim; k++) x[k] = a[k] + halt[k] * (b[k] - a[k]);
            sum1 += f(x);
            for(int k = 0; k < dim; k++) x[k] = a[k] + lat[k] * (b[k] - a[k]);
            sum2 += f(x);
        }
        double mean1 = sum1 / N, mean2 = sum2 / N, mean = (mean1 + mean2) / 2, sigma = Abs(mean1 - mean2);
        var result = (mean * V, sigma * V);
        return result;
    }
    
    public static (double, double) stratmc(Func<vector, double> f, vector a, vector b, int N, int nmin = 500){
        if(N < nmin) return plainmc(f, a, b, N);
        int dim = a.size; 
        double V = 1; 
        for(int i = 0; i < dim; i++) V *= b[i] - a[i];
        double sum = 0, sum2 = 0;
        var x = new vector(dim);
        var rnd = new Random();
        matrix sums = new matrix(2, dim), sums2 = new matrix(2, dim), vars = new matrix(2, dim);
        
        for(int i = 0; i < nmin; i++){
            for(int k = 0; k < dim; k++) x[k] = a[k] + rnd.NextDouble() * (b[k] - a[k]);
            double fx = f(x); 
            sum += fx; 
            sum2 += fx * fx;
            for(int k = 0; k < dim; k++){
                if(x[k] < (b[k] + a[k]) / 2){
                    sums[0, k] += fx; 
                    sums2[0, k] += fx * fx;
                } 
                else {
                    sums[1, k] += fx; 
                    sums2[1, k] += fx * fx;
                }
            }
        }
        
        double mean = sum / nmin, sigma = V * Sqrt(sum2 / nmin - mean * mean) / Sqrt(nmin);
        int wdim = 0; 
        double maxvar = 0;
        for(int k = 0; k < dim; k++){
            vars[0, k] = (sums2[0, k] - sums[0, k] * sums[0, k] / nmin) / nmin;
            vars[1, k] = (sums2[1, k] - sums[1, k] * sums[1, k] / nmin) / nmin;
            if(vars[0, k] > maxvar){
                maxvar = vars[0, k]; 
                wdim = k;
            }
            if(vars[1, k] > maxvar){
                maxvar = vars[1, k]; 
                wdim = k;
            }
        }
        
        vector a_new = a.copy(), b_new = b.copy();
        a_new[wdim] = (a[wdim] + b[wdim]) / 2; 
        b_new[wdim] = (a[wdim] + b[wdim]) / 2;

        int N_est = (int)Floor((N - nmin) / (1 + vars[0, wdim] / vars[1, wdim]));
        int N_up = Min(N - nmin - 2, Max(2, N_est));
        int N_down = N - nmin - N_up;

        (double int_down, double sigma_down) = stratmc(f, a, b_new, N_down, nmin);
        (double int_up, double sigma_up) = stratmc(f, a_new, b, N_up, nmin);
        double grandint = ((int_down + int_up) * (N - nmin) + V * mean * nmin) / N;
        double grandsigma = Sqrt(Pow(sigma_down * (N - nmin) / N, 2) + Pow(sigma_up * (N - nmin) / N, 2) + Pow(sigma * nmin / N, 2));
        return (grandint, grandsigma);
    }
}

