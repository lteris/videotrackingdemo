/*
 * scene_analyzer.h
 *
 *  Created on: May 28, 2009
 *      Author: liviu
 */

#ifndef SCENE_ANALYZER_H_
#define SCENE_ANALYZER_H_

#include <vector>
#include <cstdlib>
#include <ctime>
#include "digital_tree.h"

#define MAX_LIST_LEN 1000

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
			if (SameType<T, int> ) {
				memset(vals, -MAX_LIST_LEN, lines * columns * sizeof(T));
			} else if (SameType<T, bool> ) {
				memset(vals, false, lines * columns * sizeof(T));
			}
		}

		virtual ~MatchVector2D() {
			delete[] vals;
		}

		void init(int lines, int columns) {
			delete[] vals;
			vals = new T[lines * columns];
			if (SameType<T, int> ) {
				memset(vals, -MAX_LIST_LEN, lines * columns * sizeof(T));
			} else if (SameType<T, bool> ) {
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
		typename GNG_T::Node CCNode;
		typename vq::Labelizer<GNG_T>::ConnectedComponent::Nodes::iterator
				node_iter, node_end;
		mirage::img::Coordinate A, B;
		int operator()(const CCNode& model, const CCNode& graph) {
			A((int) (model.value.w[0]), (int) (model.value.w[1]));
			B((int) (graph.value.w[0]), (int) graph.value.w[1]);
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

/** CCOMPONENT class describing a connected component
 * SIMILARITY functor that computes the similarity of two connected components
 */

template<typename CCOMPONENT, typename SIMILARITY>
class RTABUMatch {
	private:

		SIMILARITY sim;

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
					this->repetitions++;
				}
		};

		class Assignment {
			private:

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
		};

		/** describes an elementary move and the associated decrease of the fitness function */
		class ElementaryMove {
			public:
				int mNode;
				int gNode;

				double fitnessDecrease;

				ElementaryMove(int mNode, int gNode) {
					this->mNode = mNode;
					this->gNode = gNode;
				}
		};

		Assignment crtConfig;
		Assignment bestConfig;
		/* latestOccupation(x, y) = last time node x in the model was assigned to node y in the graph */
		MatchVector2D<int> latestOccupation;
		/** tree that holds all the visited configurations */
		digital_tree::Tree<StoredMatch> configTree;

		CCOMPONENT* model;

		int chaotic; /* number of often repeated placements */
		int listSize;
		int stepsSinceSizeChange;
		double movingAverage;
		const int p_maxIterations;
		const double p_subOptimum; /* target similarity value */
		const int p_cycleMax;
		const int p_chaos;
		const int p_rep;
		const double p_listInc;
		const double p_listDec;

		std::vector<ElementaryMove> crtAvailableMoves;

		std::vector<mirage::img::Coordinate> modelCoordinates;
		std::vector<mirage::img::Coordinate> graphCoordinates;

	private:
		/* record time of latest occupation for the vertex to be given a new assignment */
		void makeTabu(int r) {
			latestOccupation(r, crtConfig.phi[r]) = crtConfig.time;
		}

		/** aspiration criterion - return true if by assigning m to g the similarity is better than the best one
		 * encoutered yet
		 */
		bool aspiration(const CCOMPONENT& graph, ElementaryMove& mv) {
			double pFit = crtConfig.prospectiveFitness(*model, graph, mv.mNode,
					mv.gNode);
			mv.fitnessDecrease = crtConfig.fitness(*model, graph) - pFit;
			if (pFit < bestConfig.fitness(*model, graph)) {
				return true;
			} else {
				return false;
			}
		}

		/** return true if g was assigned to m sometime during the past listSize iterations */
		bool moveIsTabu(const ElementaryMove& mv) {
			if (latestOccupation(mv.mNode, mv.gNode) >= crtConfig.time
					- listSize) {
				return true;
			} else {
				return false;
			}
		}

		/** initalize all structures when a new graph is fed to the module */
		void init(const CCOMPONENT& graph) {
			latestOccupation.init(model->nodes.size(), graph.nodes.size());
			crtConfig.init(*model, graph);
			bestConfig.init(*model, graph);
			chaotic = 0;
			listSize = 1;
			stepsSinceSizeChange = 0;
			movingAverage = 0;

			/* set the coordinates for the nodes in the graph */
			CCOMPONENT::Nodes::iterator it;
			CCOMPONENT::Node* crt;
			mirage::img::Coordinate A;
			graphCoordinates.clear();
			graphCoordinates.reserve(graph.Nodes.size());
			for (it = graph.nodes.begin(); it != graph.nodes.end(); it++) {
				graphCoordinates[(*it)->idf] = A((int) ((*it)->value.w[0]),
						(int) ((*it)->value.w[1]));
			}

			/* find an initial match */
			greedyMatch(graph, crtConfig);
			bestConfig = crtConfig;
			configTree.clear();
			configTree.resize(graph.nodes.size());
		}

		/** returns true when escape mechanism has to be performed */
		bool checkForRepetitions(const Assignment& crt) {
			stepsSinceLastChange++;
			StoredMatch *sm = configTree.findNode(crt.phi);
			if (sm != NULL) {
				int cycle_len = crt.time - sm->time;
				sm->time = crt.time;
				sm->repetitions++;
				if (sm->repetitions > p_rep) {
					chaotic++;
					if (chaotic > p_chaos) {
						chaotic = 0;
						return true;
					}
				}
				if (cycle_len < p_cycleMax) {
					movingAverage = 0.1 * cycle_len + 0.9 * movingAverage;
					listSize = (int) (listSize * p_listInc);
					stepsSinceSizeChange = 0;
				}
			} else {
				/* record the new configuration */
				configTree.insertNode(crt.phi, new StoredMatch(crt.time));
			}

			if (stepsSinceSizeChange > movingAverage) {
				listSize = (listSize * p_listDec > 1) ? listSize * p_listDec
						: 1;
				stepsSinceSizeChange = 0;
			}
			return false;
		}

		/** find best local minimum with respect to the similarity functor*/
		void greedyMatch(const CCOMPONENT& graph, Assignment& assig) {
			/* for every node in the model find the closest matching one in the graph */
			mirage::img::Coordinate M, G;

			A((int) (model.value.w[0]), (int) (model.value.w[1]));
			B((int) (graph.value.w[0]), (int) graph.value.w[1]);
			int dist = abs(A[0] - B[0]) + abs(A[1] - B[1]);

			for (int i = 0; i < modelCoordinates.size(); i++) {

			}

		}

		/** elementary moves available for the current assignement */
		void findElementaryMoves(const CCOMPONENT& graph,
				const Assignment& assig) {
			crtAvailableMoves.clear();
			/* for every node in the model s.t. phi[m]==g generate phi[m]=gi
			 * where gi is a neighbour of g in the graph
			 */
			for (int i = 0; i < crtConfig.phi.size(); i++) {
				/* get pointer to node with index phi[i] */
				CCOMPONENT::Nodes::iterator it;
				CCOMPONENT::Node* crt;
				for (it = graph.nodes.begin(); it != graph.nodes.end(); it++) {
					if ((*it)->idf == phi[i]) {
						crt = *it;
						break;
					}
				}
				int idx;
				/* find all the neighbors of crt - iterate through the edges*/
				CCOMPONENT::Edges::iterator ite;
				for (ite = crt->edges.begin(); ite != crt->edges.end(); ite++) {
					if ((*ite)->n1->idf != crt->idf) {
						idx = (*ite)->n1->idf;
					} else {
						idx = (*ite)->n2->idf;
					}
					crtAvailableMoves.push_back(ElementaryMove(i, idx));
				}
			}
		}

		/** choose best move from the crtAvailableMoves */
		void chooseBestMove(const CCOMPONENT& graph, int& m_chosen,
				int& g_chosen) {
			for (int i = 0; i < crtAvailableMoves.size(); i++) {
				/* find a move that is not tabu or satisfies the aspiration criterion */
				if (aspiration(graph, crtAvailableMoves[i]) || !moveIsTabu(
						crtAvailableMoves[i])) {
					m_chosen = crtAvailableMoves[i].mNode;
					g_chosen = crtAvailableMoves[i].gNode;
					return;
				}
			}

			/* find the best move (the ma fiteness decrease) disregarding of their tabu status */
			ElementaryMove& best = crtAvailableMoves[0];
			for (int i = 1; i < crtAvailableMoves.size(); i++) {
				if (crtAvailableMoves[i].fitnessDecrease > best.fitnessDecrease) {
					best = crtAvailableMoves[i];
				}
			}
			m_chosen = best.mNode;
			g_chosen = best.gNode;
		}

		/** escape mechanism - execute a random number of random elementary moves*/
		void doEscape(const CCOMPONENT& graph) {
			/* clean the hashing memory structure */
			configTree.clean();
			srand((unsigned) time(0));
			double random_double = (double) rand() / RAND_MAX;
			int steps = (int) (1 + (1 + random_double) * movingAverage / 2);
			/* execute steps random moves - exchange 2 units randomly */
			for (int i = 0; i < steps; i++) {
				int s1, s2;
				do {
					srand((unsigned) time(0));
					s1 = (rand() % crtConfig.phi.size());
					s2 = (rand() % crtConfig.phi.size());
				} while (s1 == s2);
				int aux1 = crtConfig.phi[s1];
				int aux2 = crtConfig.phi[s2];
				crtConfig.match(s1, aux1) = false;
				crtConfig.match(s2, aux2) = false;
				crtConfig.match(s1, aux2) = true;
				crtConfig.match(s2, aux1) = true;
				phi[s1] = aux2;
				phi[s2] = aux1;
				crtConfig.noFitness = true;
			}
			if (crtConfig.fitness(*model, graph) < bestConfig.fitness(*model,
					graph)) {
				bestConfig = crtConfig;
			}
		}

	public:
		RTABUMatch(const int maxIterations, const double subOptimum,
				const int cycleMax, const int rep, const int chaos,
				const double listInc, const double listDec) {
			this->p_maxIterations = maxIterations;
			this->p_subOptimum = subOptimum;
			this->p_cycleMax = cycleMax;
			this->p_rep = rep;
			this->p_chaos = chaos;
			this->p_listInc = listInc;
			this->p_listDec = listDec;
		}

		void setModel(CCOMPONENT* model) {
			this->model = model;
			/* set the coordinates for the nodes in the model */
			CCOMPONENT::Nodes::iterator it;
			CCOMPONENT::Node* crt;
			mirage::img::Coordinate A;
			modelCoordinates.clear();
			modelCoordinates.reserve(model->Nodes.size());
			for (it = model->nodes.begin(); it != model->nodes.end(); it++) {
				modelCoordinates[(*it)->idf] = A((int) ((*it)->value.w[0]),
						(int) ((*it)->value.w[1]));
			}
		}

		/** returns true if the similarity falls below the given threshold */
		bool operator()(const CCOMPONENT& graph) {
			/* new assignment to be performed */
			int chosen_m;
			int chosen_g;

			init(graph);

			while (current.time < p_maxIterations) {
				/* find all possible elementary moves */
				findElementaryMoves(graph, crtConfig);
				if (checkForRepetitions(crtConfig) == false) {
					/* no escape */
					chooseBestMove(chosen_m, chosen_g);
					/* record the latest occupation for the vertex in the model (chosen_m) that is to be
					 * assigned to chosen_g
					 */
					makeTabu(chosen_m);
					crtConfig.makeMove(chosen_m, chosen_g);
					crtConfig.time++;
					if (crtConfig.fitness(*model, graph) < bestConfig.fitness(
							*model, graph)) {
						bestConfig = crtConfig;
					}
				} else {
					/* execute escape mechanism */
					doEscape();
				}

				if (bestConfig.fitness(*model, graph) < p_subOptimum) {
					return true;
				}
			}

			return false;
		}
};

#endif /* SCENE_ANALYZER_H_ */
