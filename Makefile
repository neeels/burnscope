all: burnscope burnscope3 fftw3_test burnscope_fft

override CFLAGS += -Wall -O3
#override CFLAGS += -g

.PHONY: clean
clean:
	rm -f burnscope burnscope3 fftw3_test burnscope_fft

burnscope3: burnscope3.c
	$(CC) $(CFLAGS) burnscope3.c -o burnscope3 -lm -lSDL2

burnscope: burnscope.c
	$(CC) $(CFLAGS) burnscope.c -o burnscope -lm -lSDL2 -lpng

fftw3_test: fftw3_test.c
	$(CC) $(CFLAGS) fftw3_test.c -o fftw3_test -lm -lSDL2 -lfftw3

burnscope_fft: burnscope_fft.c images.h palettes.h
	$(CC) $(CFLAGS) burnscope_fft.c -o burnscope_fft -lSDL2 -lfftw3_threads -lfftw3 -lm -lpng -lsndfile

# vim: noexpandtab
