using System;
using static System.Console;
using System.Linq;
using System.Diagnostics;
using System.Threading.Tasks;

public class main{
	public class data { public int a,b; public double sum;}
	public static void harm(object obj){
		var arg = (data)obj;
		arg.sum=0;
		for(int i=arg.a;i<arg.b;i++)arg.sum+=1.0/i;
		}


	public static int Main(string[] args){
		int nthreads = 1, nterms = (int)1e8; /* default values */
		foreach(var arg in args) {
			var words = arg.Split(':');
			if(words[0]=="-threads") nthreads=int.Parse(words[1]);
			if(words[0]=="-terms"  ) nterms  =(int)float.Parse(words[1]);
			}
		WriteLine($"Terms = {nterms}, Threads = {nthreads}");

		data[] datas = new data[nthreads];
		for(int i=0;i<nthreads;i++) {
			datas[i] = new data();
			datas[i].a = 1 + nterms/nthreads*i;
			datas[i].b = 1 + nterms/nthreads*(i+1);
			WriteLine($"i = {i}, a = {datas[i].a}, b = {datas[i].b}");
			}
		datas[datas.Length-1].b=nterms+1; /* the enpoint might need adjustment */

		var threads = new System.Threading.Thread[nthreads];
		for(int i=0;i<nthreads;i++) {
			threads[i] = new System.Threading.Thread(harm); /* create a thread */
			threads[i].Start(datas[i]); /* run it with params[i] as argument to "harm" */
			}

		foreach(var thread in threads){ 
			thread.Join();
			}
		double total=0;
		foreach(var p in datas) {
			total+=p.sum;
			}
		WriteLine($"Main:total= {total}");

		var sum = new System.Threading.ThreadLocal<double>( ()=>0, trackAllValues:true);
		System.Threading.Tasks.Parallel.For( 1, nterms+1, (int i)=>sum.Value+=1.0/i );
		double totalsum=sum.Values.Sum();
		WriteLine($"Main: Sum = {totalsum}");
		return 0;
	}
}
