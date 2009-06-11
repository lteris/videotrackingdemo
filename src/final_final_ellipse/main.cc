#include "ellipse_matching.h"
#include <vq.h>
#include <stdio.h>

typedef vq::Vector<2, double> Coordinate;
typedef vq::Graph<Coordinate, int> Graph;

int main(int argc, char** argv) {

	CurvePanel cp;

	Coordinate x;

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
		cp.addPoint(x[0], x[1]);
	}
	file.close();

	double a, b, c, d, e, f;
		cp.fitStraightEllipse(a, b, c, d, e, f);
	//
	//cp.fitTiltedEllipse(a, b, c, d, e, f);

	std::vector<ControlPoint> vec;
	cp.getDrawingPoints(a, b, c, d, e, f, vec, 500);

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

	for (int i = 0; i < vec.size(); i++) {
		ofile << "1 3 0 1 0 7 50 -1 -1 0.000 1 0.0000 " << vec[i].x << " "
				<< vec[i].y << " 15 15" << " 0 0 0 0\n";
	}

	ofile.close();

	//file.close();


}

