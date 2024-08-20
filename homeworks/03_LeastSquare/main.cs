using System;
using static System.Console;
using static System.Math;

public static class main {
    public static int Main(string[] args) {
        // Read the input file name from command-line arguments
        string inputFileName = null;
        foreach (var arg in args) {
            var parts = arg.Split(':');
            if (parts[0] == "-input") {
                inputFileName = parts[1];
            }
        }
        
        if (inputFileName == null) {
            Error.WriteLine("Error: Input filename is missing.");
            return 1;
        }

        // Initialize lists to store data
        genlist<double> timeList = new genlist<double>();
        genlist<double> yValuesList = new genlist<double>();
        genlist<double> dyValuesList = new genlist<double>();

        // Read data from the input file
        using (var inputStream = new System.IO.StreamReader(inputFileName)) {
            string line;
            while ((line = inputStream.ReadLine()) != null) {
                var numbers = line.Split(' ');
                timeList.add(double.Parse(numbers[0]));
                yValuesList.add(double.Parse(numbers[1]));
                dyValuesList.add(double.Parse(numbers[2]));
            }
        }

        vector time = timeList.get_data();
        vector yValues = yValuesList.get_data();
        vector dyValues = dyValuesList.get_data();

        
        vector logY = new vector(9);
        vector relativeError = new vector(9);
        for (int i = 0; i < 9; i++) {
            logY[i] = Log(yValues[i]);
            relativeError[i] = dyValues[i] / yValues[i];
        }

        // Least squares fit
        Func<double, double>[] fittingFunctions = new Func<double, double>[] {
            x => 1,
            x => -x
        };

        
        (vector popt, matrix pcov) = fit.lsfit(fittingFunctions, time, logY, relativeError);

        // Output results
        popt.print("Optimal parameters (popt):");
        pcov.print("Covariance matrix (pcov):");
        WriteLine($"Uncertainties: {Sqrt(pcov[0, 0])} {Sqrt(pcov[1, 1])}");

        double decayConstant = Log(2) / popt[1];
        double uncertaintyDecayConstant = Log(2) / (popt[1] * popt[1]) * Sqrt(pcov[1, 1]);
        WriteLine($"starting value: a = {Exp(popt[0])}");
        WriteLine($"half life: tau = {decayConstant} ± {uncertaintyDecayConstant} days, the table value is 3.6313(14) days");

        

        return 0;
    } // Main
} // main

