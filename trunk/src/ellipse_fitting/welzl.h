/*
 * welzl.h
 *
 *  Created on: Jun 5, 2009
 *      Author: liviu
 */

#ifndef WELZL_H_
#define WELZL_H_

#include<math.h>

template<typename GRAPH>
class MinCircle {
	public:
		void operator()(const GRAPH& graph, double& xc, double& yc, double& radius) {
			/* get a point in the graph and build a circle around it */
			typename GRAPH::Nodes::iterator it;
			double x, y;

			it = graph.nodes.begin();
			xc = (*it)->value[0];
			yc = (*it)->value[1];

			it++;
			it++;
			it++;
			x = (*it)->value[0];
			y = (*it)->value[1];

			double dx = xc;
			double dy = yc;

			radius = fabs(x - dx) + fabs(y - dy);

		}
};

#endif /* WELZL_H_ */
