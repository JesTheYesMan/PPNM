using System;
using static System.Console;
using static System.Math;
using static cmath;

public class math{
	static void Main(){
		complex i = I;

	//sqrt(-1)
	complex comn = new complex(-1,0);
	Console.WriteLine("sqrt(-1):");
	Console.WriteLine($"sqrt(-1) = {cmath.sqrt(comn)}");
	Console.WriteLine($"Should equal i");
	Console.WriteLine();

	//sqrt(i)
	complex sqrti = cmath.sqrt(i);
	Console.WriteLine("sqrt(i):");
	Console.WriteLine($"sqrt(i) = {sqrti}  {sqrti}");
        Console.WriteLine($"Should equal 0.70... + 0.70... i");
	Console.WriteLine();

	//e^i
	complex ePowi = cmath.exp(i);
	Console.WriteLine("e^i:");
	Console.WriteLine($"exp(i) = {ePowi}  {ePowi}");
        Console.WriteLine($"Should equal 0.54... + 0.84... i");
	Console.WriteLine();

	//e^(i pi)
	complex ePowipi = cmath.exp(i*Math.PI);
	Console.WriteLine("e^(i pi):");
	Console.WriteLine($"exp(ipi) = {ePowipi}");
	Console.WriteLine($"Should equal -1");
	Console.WriteLine();
        

	//i^i
	complex iPowi = i.pow(i);
	Console.WriteLine("i^i:");
	Console.WriteLine($"i^i = {iPowi}");
        Console.WriteLine($"Should equal 0.20...");
	Console.WriteLine();

	//ln(i)
	complex lni = cmath.log(i);
	Console.WriteLine("ln(i):");
	Console.WriteLine($"ln(i) = {lni}");
        Console.WriteLine($"Should equal 1.57...i");
	Console.WriteLine();

	//sin(i pi)
	complex sinipi = cmath.sin(i*Math.PI);
	Console.WriteLine("sin(i pi):");
	Console.WriteLine($"sin(i pi) = {sinipi}");
	Console.WriteLine($"Should equal 11.54...i");
	Console.WriteLine();

	}
}
