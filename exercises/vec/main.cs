using System;
using static System.Math;
using static System.Console;

class main{


	static void Main(){

	// define and write the vectors
	vec a = new vec(1,2,3);
	vec b = new vec(-1,0,1);
	vec c = new vec(0.5,Math.Round(PI,2),1);
	vec d = new vec(0,0,0);
	Write("definition of the vectors for testing:\n");
	a.print("a = ");
	b.print("b = ");
	c.print("c = ");
	Write("\n");

	//testing the operations
	Write("Test of simple operations:\n");
	(2*a).print("2*a = ");
	(a+c).print("a+c = ");
	(c-b).print("c-b = ");
	(b-c).print("b-c = ");
	(-a-b).print("-a-b = ");
	Write("\n");


	//testing vector dot products
	Write($"Test of dot products(using x.dot(y) ):\n");
	Write($"a . b = {a.dot(b)}\n");
	Write($"b . a = {b.dot(a)}\n");
	Write($"a . a = {a.dot(a)}\n");
	Write($"a . d = {a.dot(d)}\n\n");


	Write($"Test of dot products(using x.dot(y) ):\n");
	Write($"a . b = {vec.dot(a,b)}\n");
        Write($"b . a = {vec.dot(b,a)}\n");
        Write($"a . a = {vec.dot(a,a)}\n");
        Write($"a . d = {vec.dot(a,d)}\n\n");


	}






}


