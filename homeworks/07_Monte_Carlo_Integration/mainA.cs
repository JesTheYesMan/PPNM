using System;
using static System.Console;
using static System.Math;

public class main{
    public static int Main(){
        Func<vector, double> circ = delegate(vector v){
            return (v.norm() < 1) ? 1 : 0;
        };

        Func<vector, double> sph = delegate(vector v){
            return (v.norm() < 1) ? Sqrt(1 - v.norm() * v.norm()) : 0;
        };

        vector a2d = new vector("-1 -1");
        vector b2d = new vector("1 1");
        vector a3d = new vector("-1 -1 -1");
        vector b3d = new vector("1 1 1");
        int numSamples = 4 * (int)1e4;
        double[] exactVals = new double[3] { PI, 2 * PI / 3, PI * PI / 4 };
        double circVal = 0, sphVal2d = 0, sphVal4d = 0;
        double circErr = 0, sphErr2d = 0, sphErr4d = 0;

        for (int i = 500; i < numSamples; i += 500){
            (circVal, circErr) = integrate.plainmc(circ, a2d, b2d, i);
            double circAbsErr = Abs(circVal - exactVals[0]);
            (sphVal2d, sphErr2d) = integrate.plainmc(sph, a2d, b2d, i);
            double sphAbsErr2d = Abs(sphVal2d - exactVals[1]);
            (sphVal4d, sphErr4d) = integrate.plainmc(sph, a3d, b3d, i);
            double sphAbsErr4d = Abs(sphVal4d - exactVals[2]);

            Error.WriteLine($"{i} {circErr} {circAbsErr} {sphErr2d} {sphAbsErr2d} {sphErr4d} {sphAbsErr4d}");
        }

        Out.WriteLine($"Circle area, N = {numSamples} exact: {exactVals[0]}, MC: {circVal}, Estimated error: {circErr}");
        Out.WriteLine($"Half sphere volume, N = {numSamples} exact: {exactVals[1]}, MC: {sphVal2d}, Estimated error: {sphErr2d}");
        Out.WriteLine($"4D half sphere volume, N = {numSamples} exact: {exactVals[2]}, MC: {sphVal4d}, Estimated error: {sphErr4d}");
        Out.WriteLine("");

        Func<vector, double> hardFn = delegate(vector v){
            return Pow(PI, -3) / (1 - Cos(v[0]) * Cos(v[1]) * Cos(v[2]));
        };

        double hardExact = 1.3932039296856768591842462603255;
        vector ha = new vector("0 0 0");
        vector hb = new vector(3); 
        for (int i = 0; i < 3; i++) hb[i] = PI;

        (double hardVal, double hardErr) = integrate.plainmc(hardFn, ha, hb, numSamples);
        int maxSamples = (int)1e6;

        Out.WriteLine($"Hard integral, N = {maxSamples}. MC: {hardVal}, est. err. {hardErr}, Exact: {hardExact}");
        return 0;
    } // Main
} // main

