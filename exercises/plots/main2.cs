using static System.Math;
class main{
public static int Main(){
for(double x=-5;x<=5.5;x+=1.0/(Pow(2,10))){
        System.Console.Write($"{x} {sfuns.gamma(x)}\n");
        }

return 0;
}//Main
}//class
