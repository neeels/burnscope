all: burnscope burnscope3

burnscope3: burnscope3.c
	gcc -g burnscope3.c -o burnscope3 -lm -lSDL

burnscope: burnscope.c
	gcc -g burnscope.c -o burnscope -lm -lSDL

# vim: noexpandtab
