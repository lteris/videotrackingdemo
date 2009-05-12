build:compile link

compile: track.h main.cc track.cc
	g++ -c -g -O3 `pkg-config --cflags gngtlib mirage-background-detection` *.cc

link: compile
	g++ -o main `pkg-config --libs gngtlib mirage bkbd-1.0` main.o track.o
