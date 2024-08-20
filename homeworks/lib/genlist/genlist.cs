using System;
using static System.Console;

public class genlist<T>{

public T[] data;
public T this[int i] {get{return data[i];} set{data[i]=value;} }
public int size = 0, capacity = 1;

public genlist(){data = new T[capacity];}
public void add (T item){
	if(size == capacity){
		T[] newdata = new T[capacity*=2];
		for(int i=0;i < size;i++)newdata[i]=data[i];
		data = newdata;
		}
	data[size] = item;
	size++;	
	}

public void remove(int i){
	if (-1 < i && i < size ){
		for(int j = i; j< size-1;j++)data[j] = data[j+1];
		size--;
		}
	else {Error.WriteLine($"Removce: Index {i} not in list of length {data.Length}");}
	}

public T[] get_data(){
	T[] _data = new T[size];
       	for(int i=0;i<size;i++)_data[i]=data[i];
	return _data;
	}
}

public class node<T>{
	public T item;
	public node<T> next;
	public node(T item){this.item=item;}
}

public class clist<T>{
        public node<T> first=null,current=null;
        public void add(T item){
                if(first==null){
                        first=new node<T>(item);
                        current=first;
                }
                else{
                        current.next = new node<T>(item);
                        current=current.next;
                }
        }
        public void start(){ current=first; }
        public void next(){ current=current.next; }
}
