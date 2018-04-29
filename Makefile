CC=g++
CFLAGS= -std=c++11 -Wextra -pedantic -g -lm -Wno-unused-variable -Wno-unused-parameter



main: main.cpp
	$(CC) $(CFLAGS) $^ -o tp2


clean:
	rm main
	rm mult
	rm *.txt

