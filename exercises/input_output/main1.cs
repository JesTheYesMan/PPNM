using static System.Math;
using System;
using static System.Console;

public class main{

	

	public static void Main(string[] args){
	foreach(var arg in args){
		var words = arg.Split(':');
		if(words[0]=="-numbers"){
			var numbers=words[1].Split(',');
			foreach(var number in numbers){
				double x = double.Parse(number);
				WriteLine($"{x}: sin({x}) = {Sin(x)}, and cos({x}) = {Cos(x)}");
				}
			}
		}
	}


}
