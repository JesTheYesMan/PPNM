using System;
using static System.Console;
using static System.Math;

public class main{
    public static int Main(){
        int n = 10, m = 100;
        System.Random random = new System.Random();
        double[] x = new double[n], y = new double[n], z = new double[m];
        for(int i = 0; i < n; i++){
            x[i] = i + 1;
            y[i] = random.NextDouble();
            Error.WriteLine($"{x[i]} {y[i]}");
        }
        for(int i = 0; i < m; i++){
            z[i] = x[0] + i * (x[n - 1] - x[0]) / (m - 1);
            Out.WriteLine($"{z[i]} {interp.linterp(x, y, z[i])} {interp.linterpInt(x, y, z[i])}");
        }
        return 0;
    }//Main
}//main
