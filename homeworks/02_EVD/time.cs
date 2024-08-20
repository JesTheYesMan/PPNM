using static System.Console;

class time{
	public static int Main(string[] args){
		int N = 1;
		foreach(string arg in args){
			var words = arg.Split(':');
			if(words[0] == "-msize")N= int.Parse(words[1]);
		}
		System.Random random = new System.Random();
		matrix A = new matrix(N);
		for(int i = 0; i < N; i++)
		for(int j = i; j < N; j++){ 
			A[i,j] = random.NextDouble();
			A[j,i] = A[i,j];
		}
		jacobi.cyclic(A);
		return 0;
	}//Main
}//time
