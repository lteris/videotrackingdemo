/*
 * rtabu.h
 *
 *  Created on: Jun 1, 2009
 *      Author: liviu
 */

#ifndef RTABU_H_
#define RTABU_H_

#include "rtabu_aux.h"

/** CCOMPONENT class describing a connected component
 * SIMILARITY functor that computes the similarity of two connected components
 */
template<typename CCOMPONENT, typename SIMILARITY>
class RTABUMatch {
	private:

		SIMILARITY sim;

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

		Assignment<CCOMPONENT, SIMILARITY> crtConfig;
		Assignment<CCOMPONENT, SIMILARITY> bestConfig;
		/* latestOccupation(x, y) = last time node x in the model was assigned to node y in the graph */
		MatchVector2D<int> latestOccupation;
		/** tree that holds all the visited configurations */
		digital_tree::Tree<StoredMatch> configTree;

		CCOMPONENT* model;

		int chaotic; /* number of often repeated placements */
		int listSize;
		int stepsSinceSizeChange;
		double movingAverage;
		int p_maxIterations;
		double p_subOptimum; /* target similarity value */
		int p_cycleMax;
		int p_chaos;
		int p_rep;
		double p_listInc;
		double p_listDec;

		std::vector<ElementaryMove> crtAvailableMoves;

		std::vector<Input> modelCoordinates;
		std::vector<Input> graphCoordinates;

	private:
		/* record time of latest occupation for the vertex to be given a new assignment */
		void makeTabu(int r) {
			latestOccupation(r, crtConfig.phi[r]) = crtConfig.time;
		}

		/** aspiration criterion - return true if by assigning m to g the similarity is better than the best one
		 * encountered yet
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

		/** find best local minimum with respect to the similarity functor*/
		void greedyMatch(const CCOMPONENT& graph, Assignment<CCOMPONENT,
				SIMILARITY>& assig) {
			/* for every node in the model find the closest matching one in the graph */
			CoordinateDistance coordDist;
			for (int i = 0; i < modelCoordinates.size(); i++) {
				int min = coordDist(modelCoordinates[i], graphCoordinates[0]);
				int minIdx = 0;
				for (int j = 1; j < graphCoordinates.size(); j++) {
					int d;
					if ((d
							= coordDist(modelCoordinates[i],
									graphCoordinates[j])) < min) {
						min = d;
						minIdx = j;
					}
				}
				assig.makeMove(i, minIdx);
			}
			std::cout << modelCoordinates.size() << " "
					<< graphCoordinates.size() << "\n";
			assig.printPhi();

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
			Input coord;

			/* set the coordinates for the nodes in the graph */
			typename CCOMPONENT::Nodes::iterator it;
			graphCoordinates.clear();
			graphCoordinates.resize(graph.nodes.size());
			for (it = graph.nodes.begin(); it != graph.nodes.end(); it++) {
				coord[0] = (int) ((*it)->value[0]);
				coord[1] = (int) ((*it)->value[1]);
				graphCoordinates[(*it)->idf] = coord;
			}

			modelCoordinates.clear();
			modelCoordinates.resize(model->nodes.size());
			for (it = model->nodes.begin(); it != model->nodes.end(); it++) {
				coord[0] = (int) ((*it)->value[0]);
				coord[1] = (int) ((*it)->value[1]);
				modelCoordinates[(*it)->idf] = coord;
			}

			for (int i = 0; i < graphCoordinates.size(); i++) {
				coord = graphCoordinates[i];
				std::cout << coord[0]<<"\n";
			}
			/* find an initial match */
			greedyMatch(graph, crtConfig);
			bestConfig = crtConfig;
			configTree.clear();
			configTree.resize(graph.nodes.size());
		}

		/** returns true when escape mechanism has to be performed */
		bool checkForRepetitions(const Assignment<CCOMPONENT, SIMILARITY>& crt) {
			stepsSinceSizeChange++;
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

		/** elementary moves available for the current assignement */
		void findElementaryMoves(const CCOMPONENT& graph, const Assignment<
				CCOMPONENT, SIMILARITY>& assig) {
			crtAvailableMoves.clear();
			/* for every node in the model s.t. phi[m]==g generate phi[m]=gi
			 * where gi is a neighbour of g in the graph
			 */
			for (int i = 0; i < crtConfig.phi.size(); i++) {
				/* get pointer to node with index phi[i] */
				typename CCOMPONENT::Nodes::iterator it;
				typename CCOMPONENT::Node* crt;
				for (it = graph.nodes.begin(); it != graph.nodes.end(); it++) {
					if ((*it)->idf == crtConfig.phi[i]) {
						crt = *it;
						break;
					}
				}
				int idx;
				/* find all the neighbors of crt - iterate through the edges*/
				typename CCOMPONENT::Edges::iterator ite;
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
			configTree.clear();
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
				crtConfig.phi[s1] = aux2;
				crtConfig.phi[s2] = aux1;
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

		virtual ~RTABUMatch() {

		}

		void setModel(CCOMPONENT* model) {
			this->model = model;
			/* set the coordinates for the nodes in the model */
			typename CCOMPONENT::Nodes::iterator it;
			typename CCOMPONENT::Node* crt;
			mirage::img::Coordinate A;
			modelCoordinates.clear();
			modelCoordinates.reserve(model->nodes.size());
			for (it = model->nodes.begin(); it != model->nodes.end(); it++) {
				modelCoordinates[(*it)->idf] = A((int) ((*it)->value[0]),
						(int) ((*it)->value[1]));
			}
		}

		/** returns true if the similarity falls below the given threshold */
		bool operator()(const CCOMPONENT& graph) {
			/* new assignment to be performed */
			int chosen_m;
			int chosen_g;

			init(graph);

			while (crtConfig.time < p_maxIterations) {
				/* find all possible elementary moves */
				findElementaryMoves(graph, crtConfig);
				if (checkForRepetitions(crtConfig) == false) {
					/* no escape */
					chooseBestMove(graph, chosen_m, chosen_g);
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
					doEscape(graph);
				}

				if (bestConfig.fitness(*model, graph) < p_subOptimum) {
					return true;
				}
			}

			return false;
		}
};

#endif /* RTABU_H_ */
