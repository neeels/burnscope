all: burnscope burnscope3 fftw3_test burnscope_fft

.PHONY: clean
clean:
	rm -f burnscope burnscope3 fftw3_test burnscope_fft

burnscope3: burnscope3.c
	gcc -Wall -O3 -g burnscope3.c -o burnscope3 -lm -lSDL2

burnscope: burnscope.c
	gcc -Wall -O3 -g burnscope.c -o burnscope -lm -lSDL2 -lpng

fftw3_test: fftw3_test.c
	gcc -Wall -O3 -g fftw3_test.c -o fftw3_test -lm -lSDL2 -lfftw3

burnscope_fft: burnscope_fft.c images.h palettes.h
	gcc -Wall -O3 -g burnscope_fft.c -o burnscope_fft -lSDL2 -lfftw3_threads -lfftw3 -lm -lpng -lsndfile

# vim: noexpandtab
