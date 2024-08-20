using System;
using static System.Console;
using static System.Math;

class main
{
    static int Main()
    {
        double t0 = 0;
        double tL = 2.3;
        double deltar = 1;
        vector r0 = new vector("-0.97000436 0.24308753 0.4662036850 0.4323657300 0 0 -0.93240737 -0.86473146 0.97000436 -0.24308753 0.4662036850 0.4323657300");

        Func<double, vector, vector> f = (t, r) =>
        {
            var w = new vector(12);
            for (int i = 0; i < 3; i++)
            {
                w[4 * i] = r[4 * i + 2];
                w[4 * i + 1] = r[4 * i + 3];
            }
            for (int i = 0; i < 3; i++)
            {
                for (int k = 0; k < 3; k++)
                {
                    if (k != i)
                    {
                        deltar = Pow(Pow(r[4 * k] - r[4 * i], 2) + Pow(r[4 * k + 1] - r[4 * i + 1], 2), -1.5);
                        w[4 * i + 2] += deltar * (r[4 * k] - r[4 * i]);
                        w[4 * i + 3] += deltar * (r[4 * k + 1] - r[4 * i + 1]);
                    }
                }
            }
            return w;
        };

        var rs = new genlist<vector>();
        ODE.driver(f, t0, r0, tL, ylist: rs);

        for (int i = 0; i < rs.size; i++)
        {
            WriteLine($"{rs[i][0]} {rs[i][1]} {rs[i][4]} {rs[i][5]} {rs[i][8]} {rs[i][9]}");
        }
        return 0;
    }
}

