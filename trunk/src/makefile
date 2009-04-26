build:compile link

compile: track.h main.cc track.cc
	g++ -c `pkg-config --cflags gngtlib mirage bkbd-1.0` *.cc

link: compile
	g++ -o main `pkg-config --libs gngtlib mirage bkbd-1.0` main.o track.o
