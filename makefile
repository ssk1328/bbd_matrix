all: run

run: 	a.out
	./a.out
	rm a.out	

a.out: 	t2.c
	gcc -o a.out bbd.c -llapack -std=c99 -fopenmp    
clean:
	 rm a.out
