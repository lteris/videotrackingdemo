/*
 * scene_analyzer.h
 *
 *  Created on: May 28, 2009
 *      Author: liviu
 */

#ifndef RTABU_AUX_H_
#define RTABU_AUX_H_

#include <vector>
#include <cstdlib>
#include <ctime>
#include <vq.h>
#include <mirage.h>
#include "digital_tree.h"

#define MAX_LIST_LEN 1000

typedef mirage::img::Coordinate Input;
typedef vq::Graph<Input, int> Graph;

template<typename T1, typename T2>
struct SameType {
		enum {
			value = false
		};
};

template<typename T>
struct SameType<T, T> {
		enum {
			value = true
		};
};

template<typename T>
class MatchVector2D {
	private:
		int lines;
		int columns;
		T* vals;
	public:
		MatchVector2D() {
			lines = columns = 0;
			vals = NULL;
		}

		MatchVector2D(int lines, int columns) {
			vals = new T[lines * columns];
			if (SameType<T, int>::value) {
				memset(vals, -MAX_LIST_LEN, lines * columns * sizeof(T));
			} else if (SameType<T, bool>::value) {
				memset(vals, false, lines * columns * sizeof(T));
			}
		}

		virtual ~MatchVector2D() {
			delete[] vals;
		}

		void init(int lines, int columns) {
			delete[] vals;
			this->lines = lines;
			this->columns = columns;
			vals = new T[lines * columns];
			if (SameType<T, int>::value) {
				memset(vals, -MAX_LIST_LEN, lines * columns * sizeof(T));
			} else if (SameType<T, bool>::value) {
				memset(vals, false, lines * columns * sizeof(T));
			}
		}

		T& operator()(int lin, int col) {
			return vals[lin * lines + col];
		}

		MatchVector2D<T>& operator=(const MatchVector2D<T>& rval) {
			if (this->lines != rval.lines || this->columns != rval.columns) {
				delete[] vals;
				this->lines = rval.lines;
				this->columns = rval.columns;
				vals = new T[lines * columns];
			}
			for (int i = 0; i < lines * columns; i++) {
				vals[i] = rval.vals[i];
			}
			return *this;
		}

};

class NodeDistance {
	public:
		typedef Graph::Node CCNode;
		mirage::img::Coordinate A, B;
		int operator()(const CCNode& model, const CCNode& graph) {
			A((int) (model.value[0]), (int) (model.value[1]));
			B((int) (graph.value[0]), (int) graph.value[1]);
			return abs(A[0] - B[0]) + abs(A[1] - B[1]);
		}
};

class CoordinateDistance {
	public:
		mirage::img::Coordinate A, B;
		int operator()(const mirage::img::Coordinate& A,
				const mirage::img::Coordinate& B) {
			return abs(A[0] - B[0]) + abs(A[1] - B[1]);
		}
};

/** compute similarity between two given graphs wrt the given match */
template<typename CCOMPONENT>
class SimFunction {
	public:
		/** match[x]=y   node x in the model corresponds to node y in the graph */
		double operator()(const CCOMPONENT& model, const CCOMPONENT& graph,
				const std::vector<int>& match) {
			//TODO
		}
};

/** the part of a match to be stored in the digital tree */
class StoredMatch {
	public:
		int time;
		int repetitions;

		StoredMatch(int time) {
			this->time = time;
			this->repetitions = 0;
		}

		StoredMatch& operator=(const StoredMatch& sm) {
			this->time = sm.time;
			repetitions++;
		}
};

template<typename CCOMPONENT, typename SIMILARITY>
class Assignment {
	private:
		SIMILARITY sim;
		double fitnessVal;
	public:
		/* match(x, y) == true node x in the model is assigned to node y in the graph */
		MatchVector2D<bool> match;
		/* phi[x] = y node x in the model is assigned to node y in the graph */
		std::vector<int> phi;
		int time;
		bool noFitness;

		Assignment() {
			noFitness = true;
		}

		double fitness(const CCOMPONENT& model, const CCOMPONENT& graph) {
			if (noFitness) {
				fitnessVal = sim(model, graph, phi);
				noFitness = false;
			}
			return fitnessVal;
		}

		void init(const CCOMPONENT& model, const CCOMPONENT& graph) {
			match.init(model.nodes.size(), graph.nodes.size());
			phi.clear();
			phi.reserve(model.nodes.size());
			time = 0;
			noFitness = true;
		}

		/** find the fitness measure for the current match with the exception of node m
		 *  in the model being
		 * assigned to node g in the graph
		 */
		double prospectiveFitness(const CCOMPONENT& model,
				const CCOMPONENT& graph, int m, int g) {
			std::vector<int> pPhi(phi);
			pPhi[m] = g;
			return sim(model, graph, pPhi);
		}

		Assignment& operator=(const Assignment& asg) {
			this->match = match;
			this->phi = asg.phi;
			this->time = asg.time;
			this->fitnessVal = asg.fitnessVal;
			this->noFitness = asg.noFitness;
			return *this;
		}

		/** assign node m in the model to node g in the graph */
		void makeMove(int m, int g) {
			match(m, phi[m]) = false;
			match(m, g) = true;
			phi[m] = g;
			noFitness = true;
		}

		void printPhi() {
			for (int i = 0; i < phi.size(); i++) {
				std::cout << phi[i] << "\n";
			}
		}
};

#endif /* SCENE_ANALYZER_H_ */
