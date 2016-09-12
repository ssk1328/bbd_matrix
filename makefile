all: run

run: 	a.out
	./a.out
	rm a.out	

a.out: 	t2.c
	gcc -o a.out bbd.c -llapack -lblas -std=c99 -fopenmp    
clean:
	 rm a.out
