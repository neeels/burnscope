all: burnscope burnscope3 fftw3_test

burnscope3: burnscope3.c
	gcc -Wall -g burnscope3.c -o burnscope3 -lm -lSDL

burnscope: burnscope.c
	gcc -Wall -g burnscope.c -o burnscope -lm -lSDL -lpng

fftw3_test: fftw3_test.c
	gcc -Wall -g fftw3_test.c -o fftw3_test -lm -lSDL -lfftw3

# vim: noexpandtab
