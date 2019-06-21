shared: src/exhaustive.c
	gcc --shared -fPIC -o sage/exhaustive.so src/exhaustive.c

examples: src/examples.c src/exhaustive.c
	gcc -o examples src/examples.c src/exhaustive.c
