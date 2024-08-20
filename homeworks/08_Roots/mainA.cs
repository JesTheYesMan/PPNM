using System;
using static System.Console;
using static System.Math;

public class main
{
    public static int Main()
    {
        // Example 1: Simple function f1 = x * exp(-x^2)
        Func<vector, vector> f1 = delegate (vector v)
        {
            vector w = new vector(1);
            w[0] = v[0] * Exp(-v[0] * v[0]);
            return w;
        };

        matrix xs = new matrix("-0.4 0.1 0.6");
        WriteLine("Testing Function: f(x) = x * exp(-x^2)");
        WriteLine("Using initial guesses: -0.4, 0.1, 0.6\n");

        for (int i = 0; i < 3; i++)
        {
            vector x = xs[i];
            root.newton root = new root.newton(f1, x);
            vector xmin = root.x;
            vector fmin = root.f;
            int steps = root.steps;
            int f_eval = root.f_eval;

            WriteLine($"Initial guess: v0 = {x[0]}");
            WriteLine($"Result: x_min = {xmin[0]}");
            WriteLine($"Function value at x_min: f(x_min) = {fmin[0]}");
            WriteLine($"Steps taken: {steps}, Function evaluations: {f_eval}\n");
        }

        // Example 2: Function with two variables f2(x, y) -> (1 + xy, 1 - x^2)
        Func<vector, vector> f2 = delegate (vector v)
        {
            vector w = new vector(2);
            w[0] = 1 + v[1] * v[0];
            w[1] = 1 - v[0] * v[0];
            return w;
        };

        matrix xs2 = new matrix("-0.6 0.01 0.01; 0 -0.01 0.6");
        WriteLine("Testing Function: (x, y) -> (1 + xy, 1 - x^2)");
        WriteLine("Using initial guesses: (-0.6, 0.01), (0, -0.01), (0.01, 0.6)\n");

        for (int i = 0; i < 3; i++)
        {
            vector x = xs2[i];
            root.newton root = new root.newton(f2, x);
            vector xmin = root.x;
            vector fmin = root.f;
            int steps = root.steps;
            int f_eval = root.f_eval;

            WriteLine($"Initial guess: v0 = ({x[0]}, {x[1]})");
            WriteLine($"Result: x_min = ({xmin[0]}, {xmin[1]})");
            WriteLine($"Function value at (x_min, y_min): f(x_min, y_min) = ({fmin[0]}, {fmin[1]})");
            WriteLine($"Steps taken: {steps}, Function evaluations: {f_eval}\n");
        }

        // Example 3: Function f3(x) = ((x-1)^2, ln(x))
        Func<vector, vector> f3 = delegate (vector v)
        {
            vector w = new vector(2);
            w[0] = Pow(v[0] - 1, 2);
            w[1] = Log(v[0]);
            return w;
        };

        matrix xs3 = new matrix("0.4 1.1 1.9");
        WriteLine("Testing Function: f(x) = ((x - 1)^2, ln(x))");
        WriteLine("Using initial guesses: 0.4, 1.1, 1.9\n");

        for (int i = 0; i < 3; i++)
        {
            vector x = xs3[i];
            root.newton root = new root.newton(f3, x);
            vector xmin = root.x;
            vector fmin = root.f;
            int steps = root.steps;
            int f_eval = root.f_eval;

            WriteLine($"Initial guess: x0 = {x[0]}");
            WriteLine($"Result: x_min = {xmin[0]}");
            WriteLine($"Function values at x_min: f(x_min) = ({fmin[0]}, {fmin[1]})");
            WriteLine($"Steps taken: {steps}, Function evaluations: {f_eval}\n");
        }

        // Example 4: Rosenbrock function
        Func<vector, double> Rf = delegate (vector v)
        {
            return Pow(1 - v[0], 2) + 100 * Pow(v[1] - v[0] * v[0], 2);
        };

        Func<vector, vector> dRf = delegate (vector v)
        {
            vector w = new vector(2);
            w[0] = -2 * (1 - v[0]) - 400 * v[0] * (v[1] - v[0] * v[0]);
            w[1] = 200 * (v[1] - v[0] * v[0]);
            return w;
        };

        matrix x0s_R = new matrix("1.8 0.2 0.2; 0.2 0.2 1.8");
        WriteLine("Testing Rosenbrock function (a = 1, b = 100)");
        WriteLine("Using initial guesses: (1.8, 0.2), (0.2, 0.2), (0.2, 1.8)\n");

        for (int i = 0; i < 3; i++)
        {
            vector x0 = x0s_R[i];
            root.newton root = new root.newton(dRf, x0);
            vector xmin = root.x;
            vector dfmin = root.f;
            int steps = root.steps;
            int f_eval = root.f_eval;
            double fmin = Rf(xmin);

            WriteLine($"Initial guess: v0 = ({x0[0]}, {x0[1]})");
            WriteLine($"Result: v_min = ({xmin[0]}, {xmin[1]})");
            WriteLine($"Function value at v_min: f(v_min) = {fmin}");
            WriteLine($"Gradient at v_min: df(v_min) = ({dfmin[0]}, {dfmin[1]})");
            WriteLine($"Steps taken: {steps}, Function evaluations: {f_eval}\n");
        }

        // Example 5: Himmelblau's function
        Func<vector, double> Hf = delegate (vector v)
        {
            return Pow(v[0] * v[0] + v[1] - 11, 2) + Pow(v[0] + v[1] * v[1] - 7, 2);
        };

        Func<vector, vector> dHf = delegate (vector v)
        {
            vector w = new vector(2);
            w[0] = 4 * v[0] * (v[0] * v[0] + v[1] - 11) + 2 * (v[0] + v[1] * v[1] - 7);
            w[1] = 2 * (v[0] * v[0] + v[1] - 11) + 2 * v[1] * (v[0] + v[1] * v[1] - 7);
            return w;
        };

        matrix x0s_H = new matrix("3.9 -2.1 -2.9 2.9; 3.9 2.1 -2.1 -2.1");
        WriteLine("Testing Himmelblau's function");
        WriteLine("Using initial guesses: (3.9, -2.1), (-2.9, 2.9), (-2.1, -2.1), (3.9, 2.1)\n");

        for (int i = 0; i < 4; i++)
        {
            vector x0 = x0s_H[i];
            root.newton root = new root.newton(dHf, x0);
            vector xmin = root.x;
            vector dfmin = root.f;
            int steps = root.steps;
            int f_eval = root.f_eval;
            double fmin = Hf(xmin);

            WriteLine($"Initial guess: v0 = ({x0[0]}, {x0[1]})");
            WriteLine($"Result: v_min = ({xmin[0]}, {xmin[1]})");
            WriteLine($"Function value at v_min: f(v_min) = {fmin}");
            WriteLine($"Gradient at v_min: df(v_min) = ({dfmin[0]}, {dfmin[1]})");
            WriteLine($"Steps taken: {steps}, Function evaluations: {f_eval}\n");
        }

        return 0;
    }
}

