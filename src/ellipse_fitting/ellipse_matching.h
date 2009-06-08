/*
 * ellipse_matching.h
 *
 *  Created on: Jun 3, 2009
 *      Author: liviu
 */

#ifndef ELLIPSE_MATCHING_H_
#define ELLIPSE_MATCHING_H_

#include "welzl.h"
#include <iostream>

/* the offset of the point from the expected track on the ellipse */
class PointDelta {
	public:
		double operator()(const double& cx, const double& cy, const double& el_a,
				const double& el_b, const int& px, const int& py) {
			double dx2 = (cx - px) * (cx - px);
			double dy2 = (cy - py) * (cy - py);
			double delta = dx2 + dy2 - el_a * el_a * el_b * el_b * (dx2 + dy2)
					/ (el_b * el_b * dx2 + el_a * el_a * dy2);
			return delta * delta;
		}
};

/* compute the value of the gradient in a certain point */
class PointGradient {
	private:
		double el_a; /* axis a of the ellipse */
		double el_b; /* axis b of the ellipse */
		double el_cx; /* x coordinate of the center */
		double el_cy; /* y coordinate of the center */

		double el_a2;
		double el_a4;
		double el_b2;
		double el_b4;
		double a2_b2;
		double dif_a2_b2;

	public:
		void updateEllipse(double a, double b, double cx, double cy) {
			this->el_a = (double) a;
			this->el_b = (double) b;
			this->el_cx = (double) cx;
			this->el_cy = (double) cy;

			el_a2 = el_a * el_a;
			el_a4 = el_a2 * el_a2;
			el_b2 = el_b * el_b;
			el_b4 = el_b2 * el_b2;
			a2_b2 = el_a2 * el_b2;
			dif_a2_b2 = el_a2 - el_b2;
		}

		/* update the ellipse parameters according to the gradient computed in point (px, py) */
		void operator()(double& ma, double& mb, double& mcx, double& mcy, int px,
				int py) {
			double dx2 = (px - el_cx) * (px - el_cx);
			double dy2 = (py - el_cy) * (py - el_cy);
			double sum_dx2dy2 = dx2 + dy2;
			double sum_dx2dy2_2 = (dx2 + dy2) * (dx2 + dy2);
			double num = el_b * el_b * dx2 + el_a * el_a * dy2;
			double num_3 = num * num * num;
			double sum_b2dx2_a2dy2 = el_b2 * dx2 + el_a2 * dy2;

			if (num != 0) {
				ma = -4 * el_a * el_b4 * dx2 * (el_b2 * dx2 + el_a2 * (-el_b2
						+ dy2)) * sum_dx2dy2_2 / num_3;
				mb = -4 * el_b * el_a4 * (el_b2 * dx2 + el_a2 * (-el_b2 + dy2))
						* sum_dx2dy2_2 * dy2 / num_3;
				mcx = 4 * (el_cx - px) * (dx2 - a2_b2 * (dx2 + dy2)
						/ sum_b2dx2_a2dy2 + dy2) * (1 - a2_b2 * dif_a2_b2 * dy2
						/ (num * num));
				mcy = 4 * (1 + a2_b2 * dif_a2_b2 * dx2 / (sum_b2dx2_a2dy2
						* sum_b2dx2_a2dy2)) * (dx2 - el_a2 * el_b2 * sum_dx2dy2
						/ sum_b2dx2_a2dy2 + dy2) * (el_cy - py);
			} else {
				ma = 0;
				mb = 0;
				mcx = 0;
				mcy = 0;
			}

		}
};

template<typename GRAPH, typename GRADIENT, typename DELTA>
class EllipseMatch {
	private:
		double step;
		double tolerance;
		int maxIter;
		DELTA pointDelta;
		GRADIENT pointGradient;
		MinCircle<GRAPH> minCircle;

		/** find the minimum radius circle that contains all the points - Weltzl algorithm*/
		void initialMatch(const GRAPH& graph, double& xc, double& yc,
				double& el_a, double& el_b) {

//			xc = 5625;
//			yc = 4050;
//			el_a = 1350;
//			el_b = 1350;
//			return;

			minCircle(graph, xc, yc, el_a);
			el_b = el_a;
		}

		/** get the cost of the current configuration */
		double cost(const GRAPH& graph, const double& xc, const double& yc,
				const double& el_a, const double& el_b) {
			/* parse the graph and add the deltas for each point */
			typename GRAPH::Nodes::iterator it;
			double ccost = 0;

			for (it = graph.nodes.begin(); it != graph.nodes.end(); it++) {
				ccost += pointDelta(xc, yc, el_a, el_b, (*it)->value[0],
						(*it)->value[1]);
			}
			return ccost;
		}

		/** compute the gradient for all the points in the graph */
		void gradient(const GRAPH& graph, double& grad_xc, double& grad_yc,
				double& grad_a, double& grad_b) {
			typename GRAPH::Nodes::iterator it;
			grad_a = 0;
			grad_b = 0;
			grad_xc = 0;
			grad_yc = 0;

			for (it = graph.nodes.begin(); it != graph.nodes.end(); it++) {
				double ma;
				double mb;
				double mcx;
				double mcy;
				pointGradient(ma, mb, mcx, mcy, (*it)->value[0],
						(*it)->value[1]);
				grad_a += ma;
				grad_b += mb;
				grad_xc += mcx;
				grad_yc += mcy;
			}

		}

	public:
		EllipseMatch(double step, double tolerance, int maxIter) {
			this->step = step;
			this->tolerance = tolerance;
			this->maxIter = maxIter;
		}

		/** find the closest ellipse to the points in the graph */
		double operator()(const GRAPH& graph, double& xc, double& yc, double& el_a,
				double& el_b) {

			double crtCost;
			initialMatch(graph, xc, yc, el_a, el_b);

			crtCost = cost(graph, xc, yc, el_a, el_b);

			for (int i = 0; i < maxIter; i++) {
				double grad_a;
				double grad_b;
				double grad_xc;
				double grad_yc;

				pointGradient.updateEllipse(el_a, el_b, xc, yc);
				gradient(graph, grad_xc, grad_yc, grad_a, grad_b);

				/* update the values for the ellipse coordinates according to the gradient */
				el_a -= step * grad_a;
				el_b -= step * grad_b;
				xc -= step * grad_xc;
				yc -= step * grad_yc;

				double ccost = cost(graph, xc, yc, el_a, el_b);
				std::cout << "cost " << ccost << "\n";
				if (fabs(ccost - crtCost) <= tolerance) {
					return ccost;
				}
				crtCost = ccost;
			}

			return crtCost;
		}

};

#endif /* ELLIPSE_MATCHING_H_ */
