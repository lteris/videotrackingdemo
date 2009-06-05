#include "ellipse_matching.h"
#include <vq.h>
#include <stdio.h>

typedef vq::Vector<2, double> Coordinate;
typedef vq::Graph<Coordinate, int> Graph;

int main(int argc, char** argv) {
	int iter;

	sscanf(argv[2], "%d", &iter);

	EllipseMatch<Graph, PointGradient, PointDelta> em(0.0000000001, 1, iter);

	Graph model;
	Coordinate x;
	double xc, yc, el_a, el_b;

	x[0] = 1;
	x[1] = 5;
	model.New(x);

	x[0] = 3;
	x[1] = 6;
	model.New(x);

	std::ifstream file;
	std::string comment;
	std::string mode;

	file.open("model.fig");
	if (!file) {
		std::cerr << "Cannot open file " << std::endl;
		::exit(1);
	}

	//	file >> std::ws;
	for (int i = 0; i < 10; i++) {
		std::getline(file, comment, '\n');
	}
	while (!file.eof()) {
		file >> x[0] >> x[1];
		model.New(x);
	}
	file.close();

	em(model, xc, yc, el_a, el_b);

	std::ofstream ofile;
	std::string fn;
	fn = argv[1];
	fn += "1";
	ofile.open(fn.c_str());

	file.open(argv[1]);
	while (!file.eof()) {
		std::getline(file, comment, '\n');
		ofile << comment << "\n";
	}
	file.close();
	ofile << "1 3 0 1 0 7 50 -1 -1 0.000 1 0.0000 ";
	ofile << (int) xc << " " << (int) yc << " " << (int) el_a << " "
			<< (int) el_b << " 0 0 0 0\n";

	file.close();
	ofile.close();
}

