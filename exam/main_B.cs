using static System.Console;
using static System.Math;
using System;


public class main_B{
    
    	static (matrix N, matrix H) SetupMatrices(vector alpha)
    	{
        int n = alpha.size;
        int matrix_size=n;
        matrix H = new matrix(matrix_size);
        matrix N = new matrix(matrix_size);

        for (int i = 0; i < n; i++)
        {
            for (int j = 0; j < n; j++)
            {
                N[i, j] = 1.0/4.0 *Sqrt(PI)*Pow((alpha[i]+alpha[j]),-3.0/2.0);
                H[i, j] = -3.0/2.0 *Pow(Sqrt(PI*alpha[i]*alpha[j]*(alpha[i]+alpha[j])),-5.0/2.0) + 1.0/2.0*Pow((alpha[i]+alpha[j]),-1.0);
            }
        }
//    	N.print();
//    	H.print();
    	return (N,H);
    	}
    
    
    
	public static double ObjectiveFunction(vector alpha)
	{
    	(matrix N, matrix H) = SetupMatrices(alpha);
    	var (eigenvalues, _) = main_A.generalised_eigenvalues_solver(H, N);
    	return eigenvalues[0];
	}
    
    

public static vector MinimizeAlpha()
{
    // Start with an initial guess for the alphas
    vector initialGuess = new vector(1.0,2.0,3.0,4.0); // Example initial guess
    
    // Use the qnewton optimizer to minimize the objective function
    min.qnewton optimizer = new min.qnewton(ObjectiveFunction, initialGuess, acc: 0.001);
    
    // Extract the optimized alphas and the corresponding minimum energy
    vector optimalAlphas = optimizer.x;
    double optimalEnergy = optimizer.f;
    
    Console.WriteLine("Optimal Alphas: " + optimalAlphas);
    Console.WriteLine("Optimal Ground State Energy: " + optimalEnergy);
    
    return optimalAlphas;
}



    public static void Main()
    {
        
        
        vector alpha = MinimizeAlpha();  // Example choice of variational parameters
        (matrix H, matrix N) = SetupMatrices(alpha);
        
        (vector epsilon, matrix c)=main_A.generalised_eigenvalues_solver(H, N);
        
        WriteLine("Ground State Energy: " + epsilon[0]);
        WriteLine("Ground State Vector: " + c[0]);
    }
}


