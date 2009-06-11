build:compile link

compile: main.cc
	g++ -c -O3 `pkg-config --cflags mirage vq gngtlib` *.cc

link: compile
	g++ -o main `pkg-config --libs vq mirage gngtlib` main.o
