using System;
using static System.Console;
using static System.Math;

namespace root
{
    public class newton
    {
        public static readonly double ε = Pow(2, -26), λmin = Pow(2, -10);

        public readonly Func<vector, vector> F;
        public vector x, f;
        public int steps, f_eval;
        public bool status;

        // Constructor
        public newton(Func<vector, vector> func, vector x0, int max_steps = 9999, double acc = 1e-4)
        {
            F = func;
            steps = 0;
            f_eval = 1;
            x = x0.copy();
            f = F(x0);
            status = false;
            vector Dx = null;

            do
            {
                steps++;
                matrix J = jacobian(x, f);
                Dx = QRGS.solve(J, -f);
                double lambda = 1;
                vector f1 = F(x + Dx);
                f_eval++;

                while (f1.norm() > (1 - lambda / 2) * f.norm() && λmin < lambda)
                {
                    lambda /= 2;
                    f1 = F(x + lambda * Dx);
                    f_eval++;
                }

                x += lambda * Dx;
                f = f1;

            } while (f.norm() >= acc && Dx.norm() >= ε * x.norm() && steps < max_steps);

            if (f.norm() < acc) status = true;
        }

        matrix jacobian(vector x, vector f0 = null)
        {
            if (f0 == null)
            {
                f0 = F(x);
                f_eval++;
            }

            int m = x.size, n = f0.size;
            matrix jac = new matrix(n, m);
            vector x2 = x.copy();

            for (int i = 0; i < m; i++)
            {
                double dx = ε * Max(1, Abs(x2[i]));
                x2[i] += dx;
                vector df = F(x2) - f0;
                f_eval++;
                jac[i] = df / dx;
                x2[i] = x[i];
            }

            return jac;
        }
    }
}

