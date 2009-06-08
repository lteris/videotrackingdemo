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
		double operator()(const double cm, const double cn, const double el_a,
				const double el_b, const int px, const int py) {
			double error = el_b * el_b * (cm - px) * (cm - px) + el_a * el_a
					* ((cn - py) * (cn - py) - el_b * el_b);
			return fabs(error);
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
			double mx2 = (el_cx - px) * (el_cx - px);
			double ny2 = (el_cy - py) * (el_cy - py);
			double num = el_b * el_b * mx2 + el_a * el_a * (ny2 - el_b * el_b);
			double num_abs = fabs(num);

			ma = 2 * el_a * num * (-el_b * el_b + ny2) / num_abs;
			mb = 2 * el_b * num * (-el_a * el_a + mx2) / num_abs;
			mcx = 2 * el_b * el_b * (el_cx - px) * num / num_abs;
			mcy = 2 * el_a * el_a * (el_cy - py) * num / num_abs;
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

			//			5250 3975 2112 2112

			xc = 5250;
			yc = 3975;
			el_a = 20;//2112;
			el_b = 20;//2112;
			//			minCircle(graph, xc, yc, el_a);
			//			el_b = el_a;
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
				if (ccost > crtCost) {
					return ccost;
				}
				crtCost = ccost;
			}

			return crtCost;
		}

};

#endif /* ELLIPSE_MATCHING_H_ */
