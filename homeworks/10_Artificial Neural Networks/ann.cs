using System;
using static System.Math;

namespace ann {

    public class interpol {
        public readonly int n; /* number of neurons */
        public bool train_status;
        Func<double, double> f = x => x * Exp(-x * x); /* activation function */
        Func<double, double> df = x => (1 - 2 * x * x) * Exp(-x * x); /* derivative of function */
        Func<double, double> If = x => -Exp(-x * x) / 2; /* antiderivative of function */
        public vector p; /* network parameters */
        public double val; /* value of cost function after training */

        // constructors
        public interpol(int m) {
            n = m;
            p = new vector(3 * n);
            for (int i = 0; i < n; i++) {
                p[3 * i] = 0;
                p[3 * i + 1] = 1;
                p[3 * i + 2] = 1;
            }
        }

        // Response
        double response(double x, vector v) {
            double res = 0;
            for (int i = 0; i < n; i++) res += v[3 * i + 2] * f((x - v[3 * i]) / v[3 * i + 1]);
            return res;
        }

        public double response(double x) {
            return response(x, p);
        }

        public double dresponse(double x) {
            double res = 0;
            for (int i = 0; i < n; i++) res += p[3 * i + 2] * df((x - p[3 * i]) / p[3 * i + 1]) / p[3 * i + 1];
            return res;
        } // dresponse

        public double Iresponse(double a, double b) {
            double res = 0;
            for (int i = 0; i < n; i++) res += p[3 * i + 2] * p[3 * i + 1] * (If((b - p[3 * i]) / p[3 * i + 1]) - If((a - p[3 * i]) / p[3 * i + 1]));
            return res;
        } // Iresponse

        // train
        public void train_interp(vector x, vector y, double acc = 1e-10, vector pmin = null, vector pmax = null) {
            Func<vector, double> cost = delegate (vector v) {
                double cst = 0;
                for (int i = 0; i < x.size; i++) cst += Pow(response(x[i], v) - y[i], 2);
                return cst;
            };

            if (pmin != null && pmax != null) { /* for global optimization */
                if (pmin.size != 3 * n) throw new ArgumentException($"train_interp: invalid size of parameter space: {pmin.size} != 3*dim = {3 * n}.");
                if (pmax.size != 3 * n) throw new ArgumentException($"train_interp: invalid size of parameter space: {pmax.size} != 3*dim = {3 * n}.");
                min.glo_sto minip = new min.glo_sto(cost, pmin, pmax, 5, acc: acc, d: 5);
                train_status = minip.status;
                p = minip.x;
                val = minip.f;
            } else if (pmin == null && pmax == null) { /* purely local optimization */
                min.downhill_sim minip = new min.downhill_sim(cost, p, acc: acc);
                train_status = minip.status;
                p = minip.x;
                val = minip.f;
            } else throw new ArgumentException($"train_interp: Must provide both pmin and pmax or neither.");
        }
    } // ann

    public class diff_eq {
        public readonly int n; /* number of hidden nodes */
        public bool train_status;
        Func<double, double> f = x => x * Exp(-x * x); /* activation function */
        Func<double, double> df = x => (1 - 2 * x * x) * Exp(-x * x); /* derivative of function */
        Func<double, double> ddf = x => -2 * x * (3 - 2 * x * x) * Exp(-x * x); /* double derivative */
        Func<double, double> If = x => -Exp(-x * x) / 2; /* antiderivative of function */
        public vector p; /* network parameters */

        // constructors
        public diff_eq(int m) {
            n = m;
            p = new vector(3 * n);
            for (int i = 0; i < n; i++) {
                p[3 * i] = 0;
                p[3 * i + 1] = 1;
                p[3 * i + 2] = 1;
            }
        }

        // Response
        double response(double x, vector v) {
            double res = 0;
            for (int i = 0; i < n; i++) res += v[3 * i + 2] * f((x - v[3 * i]) / v[3 * i + 1]);
            return res;
        }

        public double response(double x) {
            return response(x, p);
        }

        public double dresponse(double x) {
            return dresponse(x, p);
        }

        double dresponse(double x, vector v) {
            double res = 0;
            for (int i = 0; i < n; i++) res += v[3 * i + 2] * df((x - v[3 * i]) / v[3 * i + 1]) / v[3 * i + 1];
            return res;
        } // dresponse

        public double ddresponse(double x) {
            return ddresponse(x, p);
        }

        double ddresponse(double x, vector v) {
            double res = 0;
            for (int i = 0; i < n; i++) res += v[3 * i + 2] * ddf((x - v[3 * i]) / v[3 * i + 1]) / Pow(v[3 * i + 1], 2);
            return res;
        } // ddresponse

        public double Iresponse(double a, double b) {
            double res = 0;
            for (int i = 0; i < n; i++) res += p[3 * i + 2] * p[3 * i + 1] * (If((b - p[3 * i]) / p[3 * i + 1]) - If((a - p[3 * i]) / p[3 * i + 1]));
            return res;
        } // Iresponse

        public void train(
            Func<double[], double, double> Phi, /* diff equation = 0 */
            double a, double b, /* interval */
            double c, /* starting point */
            double[] Fc, /* initial values */
            double alp = 100,
            double bet = 100,
            double acc = 1e-10) {

            Func<vector, double> cost = delegate (vector I) {
                double res = alp * Pow(Fc[0] - response(c, I), 2) + bet * Pow(Fc[1] - dresponse(c, I), 2);
                Func<double, double[]> y = delegate (double t) {
                    vector vres = new vector(3);
                    vres[0] = response(t, I);
                    vres[1] = dresponse(t, I);
                    vres[2] = ddresponse(t, I);
                    return vres;
                };
                Func<double, double> Phi2 = x => Pow(Phi(y(x), x), 2);
                res += integrate.adint(Phi2, a, b).Item1;
                return res;
            };

            min.downhill_sim minip = new min.downhill_sim(cost, p, acc: acc);
            train_status = minip.status;
            p = minip.x;
        } // train

    } // diff_nn

} // namespace ann

