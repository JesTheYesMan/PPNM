using System;
using static System.Console;
using static System.Math;

public class interp {
    // Linear interpolation function
    public static double linterp(double[] x, double[] y, double z, int i = -1) {
        if (i == -1) {
            i = binsearch(x, z); // Perform binary search if index not provided
        }
        double dx = x[i + 1] - x[i];
        if (!(dx > 0)) throw new Exception("linterp: x-array is not properly ordered.");
        double dy = y[i + 1] - y[i];
        return y[i] + dy / dx * (z - x[i]);
    } // linterp

    // Binary search for index where z fits
    public static int binsearch(double[] x, double z) {
        if (!(x[0] <= z && x[x.Length - 1] >= z)) {
            throw new ArgumentException($"binsearch: z = {z} is outside the interval [{x[0]}, {x[x.Length - 1]}].");
        }
        int i = 0, j = x.Length - 1; // Endpoints
        while (j - i > 1) {
            int mid = (i + j) / 2;
            if (z > x[mid]) i = mid; else j = mid;
        }
        return i;
    } // binsearch

    // Integral of the linear interpolation
    public static double linterpInt(double[] x, double[] y, double z) {
        int i = binsearch(x, z);
        double res = 0;
        for (int j = 0; j < i; j++) {
            res += (x[j + 1] - x[j]) * (y[j + 1] + y[j]) / 2; // Add integrals of each trapezoid
        }
        res += (z - x[i]) * (linterp(x, y, z, i) + y[i]) / 2;
        return res;
    } // linterpInt

    public class qspline {
        public vector x, y, b, c;

        public qspline(vector xs, vector ys) {
            if (xs.size != ys.size) {
                throw new ArgumentException($"qspline: Data sizes incompatible. x: {xs.size}, y: {ys.size}");
            }
            x = xs.copy();
            y = ys.copy();
            int n = x.size;
            vector p = new vector(n - 1), dx = new vector(n - 1);
            b = new vector(n - 1);
            c = new vector(n - 1);

            for (int i = 0; i < n - 1; i++) {
                dx[i] = x[i + 1] - x[i];
                p[i] = (y[i + 1] - y[i]) / dx[i];
            }
            for (int i = 0; i < n - 2; i++) {
                c[i + 1] = (p[i + 1] - p[i] - c[i] * dx[i]) / dx[i + 1];
            }
            c[n - 2] /= 2;
            for (int j = n - 3; j >= 0; j--) {
                c[j] = (p[j + 1] - p[j] - c[j + 1] * dx[j + 1]) / dx[j];
            }
            for (int i = 0; i < n - 1; i++) {
                b[i] = p[i] - c[i] * dx[i];
            }
        } // constructor

        public double evaluate(double z) {
            int i = binsearch(x, z);
            double res = y[i] + b[i] * (z - x[i]) + c[i] * Pow(z - x[i], 2);
            return res;
        } // evaluate

        public double derivative(double z) {
            int i = binsearch(x, z);
            double res = b[i] + 2 * c[i] * (z - x[i]);
            return res;
        } // derivative

        public double integral(double z) {
            int i = binsearch(x, z);
            double res = 0;
            for (int j = 0; j < i; j++) {
                double T12 = (x[j + 1] - x[j]) * y[j] + b[j] * Pow(x[j + 1] - x[j], 2) / 2;
                double T34 = c[j] * Pow(x[j + 1] - x[j], 3) / 3;
                res += T12 + T34;
            }
            res += (z - x[i]) * y[i] + b[i] * Pow(z - x[i], 2) / 2 + c[i] * Pow(z - x[i], 3) / 3;
            return res;
        } // integral
    } // qspline

    public class cspline {
        public vector x, y, b, c, d;

        public cspline(vector xs, vector ys) {
            if (xs.size != ys.size) {
                throw new ArgumentException($"cspline: Data sizes incompatible. x: {xs.size}, y: {ys.size}");
            }
            x = xs.copy();
            y = ys.copy();
            int n = x.size;

            // Construct all relevant vectors and matrices
            vector h = new vector(n - 1);
            vector p = new vector(n - 1);
            vector B = new vector(n);
            b = new vector(n);
            c = new vector(n - 1);
            d = new vector(n - 1);
            matrix A = new matrix(n);

            for (int i = 0; i < n - 1; i++) {
                h[i] = x[i + 1] - x[i];
                p[i] = (y[i + 1] - y[i]) / h[i];
            }

            // Set up boundary conditions
            A[0, 0] = 2;
            A[n - 1, n - 1] = 2;
            A[0, 1] = 1;
            B[0] = 3 * p[0];
            B[n - 1] = 3 * p[n - 2];
            A[n - 1, n - 2] = 1;

            // Fill matrix A and vector B for solving
            for (int i = 0; i < n - 2; i++) {
                A[i + 1, i + 1] = 2 * h[i] / h[i + 1] + 2;
                A[i + 1, i + 2] = h[i] / h[i + 1];
                A[i + 1, i] = 1;
                B[i + 1] = 3 * (p[i] + p[i + 1] * h[i] / h[i + 1]);
            }

            // Solve tridiagonal system for b
            b = linsol.TriDiagSol(A, B);

            // Compute c and d
            c[0] = 0;
            for (int i = 0; i < n - 2; i++) {
                c[i + 1] = (-2 * b[i + 1] - b[i + 2] + 3 * p[i + 1]) / h[i + 1];
                d[i] = (b[i] + b[i + 1] - 2 * p[i]) / (h[i] * h[i]);
            }
            d[n - 2] = -c[n - 2] / (3 * h[n - 2]);
        } // constructor

        public double evaluate(double z) {
            int i = binsearch(x, z);
            double res = y[i] + b[i] * (z - x[i]) + c[i] * Pow(z - x[i], 2) + d[i] * Pow(z - x[i], 3);
            return res;
        } // evaluate

        public double derivative(double z) {
            int i = binsearch(x, z);
            double res = b[i] + 2 * c[i] * (z - x[i]) + 3 * d[i] * Pow(z - x[i], 2);
            return res;
        } // derivative

        public double integral(double z) {
            int i = binsearch(x, z);
            double res = 0;
            for (int j = 0; j < i; j++) {
                double T12 = (x[j + 1] - x[j]) * y[j] + b[j] * Pow(x[j + 1] - x[j], 2) / 2;
                double T34 = c[j] * Pow(x[j + 1] - x[j], 3) / 3 + d[j] * Pow(x[j + 1] - x[j], 4) / 4;
                res += T12 + T34;
            }
            res += (z - x[i]) * y[i] + b[i] * Pow(z - x[i], 2) / 2 + c[i] * Pow(z - x[i], 3) / 3 + d[i] * Pow(z - x[i], 4) / 4;
            return res;
        } // integral
    } // cspline

    public class vcspline { // Vector of cubic splines (multidimensional)
        public vector x;
        public vector[] y;
        public cspline[] splines;
        public int numPoints, dimension;

        public vcspline(vector xs, vector[] ys) {
            if (ys[0].size != xs.size) {
                throw new ArgumentException($"vcspline: x ({xs.size}) and y ({ys[0].size}) must have the same size.");
            }
            numPoints = xs.size;
            dimension = ys.Length;
            x = xs.copy();
            y = new vector[dimension];
            splines = new cspline[dimension];
            for (int i = 0; i < dimension; i++) {
                y[i] = ys[i].copy();
                splines[i] = new cspline(x, y[i]);
            }
        } // constructor

        public vector evaluate(double z) {
            vector result = new vector(numPoints);
            for (int i = 0; i < dimension; i++) {
                result[i] = splines[i].evaluate(z);
            }
            return result;
        } // evaluate

        public vector derivative(double z) {
            vector result = new vector(numPoints);
            for (int i = 0; i < dimension; i++) {
                result[i] = splines[i].derivative(z);
            }
            return result;
        } // derivative

        public vector integral(double z) {
            vector result = new vector(numPoints);
            for (int i = 0; i < dimension; i++) {
                result[i] = splines[i].integral(z);
            }
            return result;
        } // integral
    } // vcspline
} // interp

