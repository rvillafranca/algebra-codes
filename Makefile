shared: src/exhaustive.c
	gcc --shared -fPIC -Ofast -o sage/exhaustive.so src/exhaustive.c

examples: src/examples.c src/exhaustive.c
	gcc -o examples src/examples.c src/exhaustive.c
