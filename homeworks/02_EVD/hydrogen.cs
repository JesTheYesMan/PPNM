using System;
using static System.Console;
using static System.Math;

public static class hydrogen {
    public static int Main(string[] args) {
        // Default parameters
        double rmax = 10;
        double dr = 0.1;
        int EnergyLevel = 0;
        string Wavef_conv = "conv";
        
        // Parse command-line arguments
        foreach (string arg in args) {
            var parts = arg.Split(':');
            if (parts[0] == "-rmax") {
                rmax = double.Parse(parts[1]);
            } else if (parts[0] == "-dr") {
                dr = double.Parse(parts[1]);
            } else if (parts[0] == "-EnergyLevel") {
                EnergyLevel = int.Parse(parts[1]);
            } else if (parts[0] == "-Wavef_conv") {
                Wavef_conv = parts[1];
            }
        }

        int N = (int)Floor(rmax / dr) - 1;

        // Initialize vectors
        vector K_diag = new vector(N);
        vector K_offdiag = new vector(N - 1);
        vector rs = new vector(N);
        vector Vs = new vector(N);

        // Construct matrices
        for (int i = 0; i < N; i++) {
            K_diag[i] = -2;
            rs[i] = (i + 1) * dr;
            Vs[i] = -1.0 / rs[i];
            if (i != N - 1) {
                K_offdiag[i] = 1;
            }
        }

        matrix K = -0.5 * (
            matrix.diag(K_diag) +
            matrix.diag(K_offdiag, 1) +
            matrix.diag(K_offdiag, -1)
        ) * Pow(dr, -2);
        
        matrix V = matrix.diag(Vs);
        matrix H = K + V; 

        // Compute eigenvalues and eigenvectors
        (vector Energies, matrix EigenVectors) = jacobi.cyclic_opt(H);

        // Ensure positive eigenvector components
        for (int i = 0; i < 3; i++) {
            if (EigenVectors[i][1] < 0) {
                EigenVectors[i] *= -1;
            }
        }

        // Output results based on format
        if (Wavef_conv == "conv") {
            WriteLine($"{rmax} {dr} {Energies[EnergyLevel]}");
        } else if (Wavef_conv == "Wavef") {
            for (int i = 0; i < N; i++) {
                WriteLine($"{rs[i]} {EigenVectors[0][i] / Sqrt(dr)} {EigenVectors[1][i] / Sqrt(dr)} {EigenVectors[2][i] / Sqrt(dr)}");
            }
        }

        return 0;
    }
}
