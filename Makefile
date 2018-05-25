CC=g++
CFLAGS= -std=c++11 -Wextra -Wall -pedantic -g -lm -Wno-unused-variable -Wno-unused-parameter



main: main.cpp util.cpp ppmloader.cpp mediciones.cpp clasificador.cpp
	$(CC) $(CFLAGS) $^ -o tp2


clean:
	rm main
	rm *.txt

