all: run

run: 	a.out
	./a.out
	rm a.out	

a.out: 	t2.c
	gcc -o a.out bbd.c -llapack -lblas -std=c99 -fopenmp    

gen_data: 
	gcc -o a.out gen_test.c -llapack -std=c99	
	./a.out

read_data: 	gen_data
	gcc -o a.out read_test.c bbdf.c -llapack -lblas -std=c99 -fopenmp 
	./a.out

clean:
	 rm a.out
