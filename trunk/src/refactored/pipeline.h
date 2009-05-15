/*   This file is part of pipeliner lib
 *
 *   Copyright (C) 2009,  Supelec
 *
 *   Author : Liviu Teris
 *
 *   Contributor :
 *
 *   This library is free software; you can redistribute it and/or
 *   modify it under the terms of the GNU Lesser General Public
 *   License as published by the Free Software Foundation; either
 *   version 2.1 of the License, or any later version.
 *   
 *   This library is distributed in the hope that it will be useful,
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 *   Lesser General Public License for more details.
 *   
 *   You should have received a copy of the GNU Lesser General Public
 *   License along with this library; if not, write to the Free Software
 *   Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 *
 *   Contact : Herve.Frezza-Buet@supelec.fr
 *
 */


#ifndef PIPELINE_H_
#define PIPELINE_H_

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <string>
#include <cc++/thread.h>
#include <queue>

namespace pipeline {
#define MAX_FRAMES_IN_PIPE 100000	/* the sequence number counter resets at this number */

	/**
	 * @short container for images passed between frames
	 * CONTENT_TYPE needs to implement the = operator
	 */
	template<typename CONTENT_TYPE>
	class Queueable {
		public:
			CONTENT_TYPE data;
			int seqNo;

		public:
			Queueable() {
				this->seqNo = 0;
			}

			Queueable(const CONTENT_TYPE& data, int seqNo) {
				this->seqNo = seqNo;
				this->data = data;
			}

			Queueable& operator=(const Queueable& cont) {
				this->seqNo = cont.seqNo;
				this->data = cont.data;
			}
	};

	/**
	 * @short position of stage in pipe
	 */
	enum {
		FIRST, LAST, MIDDLE
	};

	/**
	 * @short dummy type to generate 2 types based on type parameters
	 */
	template<bool flag> class DummyType {
			enum {
				val = flag
			};
	};

	/**
	 * @short type comparison classes - true class
	 */
	struct TrueType {
			enum {
				value = true
			};
	};

	/**
	 * @short type comparison classes - false class
	 */
	struct FalseType {
			enum {
				value = false
			};
	};

	/**
	 * @short type comparison class
	 */
	template<typename T1, typename T2>
	struct SameType {
			typedef FalseType value;
	};

	template<typename T>
	struct SameType<T, T> {
			typedef TrueType value;
	};

	/**
	 * @short empty type list
	 */
	struct NullList {
			enum {
				length = 0
			};
	};

	/**
	 * @short all the information needed for a stage in the pipeline: the type of the output, the class that implements the () operator
	 * as required by the GeneralStage class and the number of the stage that holds the writing end of the queue that goes in the 
         * current stage (if any)
	 * OUT_FRAME must implement operator=(const OUT_FRAME&)
	 */
	template<typename OUT_FRAME, typename COMPUTATION, int Q_FROM = 0>
	struct Types4Stage {
			typedef OUT_FRAME out_frame_type;
			typedef COMPUTATION computation;
			enum {
				q_from = Q_FROM
			};
	};

	/**
	 * @short generic interface for a stage in the pipeline:
	 * Q_TYPE = type for the elements stored in the queue
	 * COMPUTE = class that overloads the () operator for the parameters
	 * compute(IN_FRAME*, OUT_FRAME*&, Q_TYPE&) if a bypass queue is involved or
	 * compute(IN_FRAME*, OUT_FRAME*&) if there is no incoming queue for that stage
	 */
	template<typename IN_FRAME, typename OUT_FRAME, typename NEXT_STAGE,
			typename Q_TYPE, typename COMPUTATION, int POSITION_IN_PIPE =
					MIDDLE>
	class GenericStage: public ost::Thread {
		public:
			int crtSeqNo; /* seq # corrensponding to the data in the front buffer */
			int bbSeqNo; /* seq # assoc. with the content of the back buffer     */
			IN_FRAME* frontBuffer;
			IN_FRAME* backBuffer;
			OUT_FRAME* outFrame; /* the output of the compute operation */

			COMPUTATION compute; /* functor */

			std::queue<Queueable<Q_TYPE> >* readBpQueue; /*bypass queue the stage reads from (if any) */
			std::queue<Queueable<OUT_FRAME> >* writeBpQueue; /*bypass queue the stage writes to (if any) */

			/* synchronization vars */
			bool wait4Input; /* thread blocked waiting on input */
			bool freshInput; /* set to true when the back buffer is refreshed */
			ost::Mutex bbMutex; /* guard access to the back buffer */

			pthread_mutex_t condMutex; /* mutex to be used with the conditional variable */
			pthread_cond_t syncCond; /* conditional variable */

		protected:
			volatile bool stopFlag;
			NEXT_STAGE* next;
			std::string debug;

		private:
			void initBuffers() {
				switch (POSITION_IN_PIPE) {
				case LAST:
					outFrame = NULL;
					frontBuffer = new IN_FRAME();
					backBuffer = new IN_FRAME();
					break;
				case FIRST:
					frontBuffer = backBuffer = NULL;
					outFrame = new OUT_FRAME();
					break;
				case MIDDLE:
					/* fallthru */
				default:
					frontBuffer = new IN_FRAME();
					backBuffer = new IN_FRAME();
					outFrame = new OUT_FRAME();
				};
			}

			void write2Next(DummyType<true> ) {

				next->writeBackBuffer(outFrame, crtSeqNo);
				pthread_mutex_lock(&(next->condMutex));
				if (next->wait4Input) { /* the next thread is waiting for input */
					next->wait4Input = false;
					next->freshInput = false;
					pthread_cond_signal(&(next->syncCond));
				} else { /* the next thread is computing */
					next->freshInput = true;
				}
				pthread_mutex_unlock(&(next->condMutex));
			}

			void write2Next(DummyType<false> ) {
				/* function does nothing for last stage */
			}

			/* call functor when there's a queue involved */
			void callCompute(DummyType<false> ) {
				compute(frontBuffer, outFrame, (readBpQueue->front()).data);
			}

			/* call functor without a queue */
			void callCompute(DummyType<true> ) {
				compute(frontBuffer, outFrame);
			}

		public:
			GenericStage() {
				initBuffers();
				wait4Input = true;
				freshInput = false;
				stopFlag = false;

				crtSeqNo = 0;

				this->next = NULL;

				pthread_mutex_init(&condMutex, NULL);
				pthread_cond_init(&syncCond, NULL);

				/* create readQ if there is a need for one */
				if (SameType<Q_TYPE, NullList>::value::value) {
					readBpQueue = NULL;
				} else {
					readBpQueue = new std::queue<Queueable<Q_TYPE> >();
				}
				writeBpQueue = NULL;

				this->debug = "poster";
			}

			GenericStage(std::string debug) {
				initBuffers();
				wait4Input = true;
				freshInput = false;
				stopFlag = false;

				crtSeqNo = 0;

				this->next = NULL;

				pthread_mutex_init(&condMutex, NULL);
				pthread_cond_init(&syncCond, NULL);

				/* create readQ if there is a need for one */
				if (SameType<Q_TYPE, NullList>::value::value) {
					readBpQueue = NULL;
				} else {
					readBpQueue = new std::queue<Queueable<Q_TYPE> >();
				}
				writeBpQueue = NULL;

				this->debug = debug;
			}

			virtual ~GenericStage() {
				delete frontBuffer;
				delete backBuffer;
				delete outFrame;
				delete readBpQueue;

				pthread_cond_destroy(&syncCond);
				if (pthread_mutex_destroy(&condMutex) != 0) {
					pthread_mutex_unlock(&condMutex);
					pthread_mutex_destroy(&condMutex);
				}

			}

			void setNext(NEXT_STAGE* next) {
				this->next = next;
			}

			void setWriteQ(std::queue<Queueable<OUT_FRAME> >* writeBpQueue) {
				this->writeBpQueue = writeBpQueue;
			}

			std::queue<Queueable<Q_TYPE> >* getReadQ() {
				return this->readBpQueue;
			}

			void run() {
				while (true) {
					waitNewInput();
					if (stopFlag) {
						break;
					}
					swapBuffers();
					flushQueue();

					callCompute(DummyType<
							SameType<Q_TYPE, NullList>::value::value> ());

					if (POSITION_IN_PIPE == FIRST) {
						/* generate a new crt number */
						crtSeqNo = (crtSeqNo++) % MAX_FRAMES_IN_PIPE;
					}

					if (writeBpQueue != NULL) {
						/* sends the output to the bypass stage - no sync required when writing */
						Queueable<OUT_FRAME> elem(*outFrame, crtSeqNo);
						writeBpQueue->push(elem);
					}
					write2Next();
				}
			}

			void stop() {
				stopFlag = true;
				if (POSITION_IN_PIPE != FIRST) {
					/* first stage is never blocked waiting for the previous stage */
					pthread_mutex_lock(&condMutex);
					if (wait4Input) { /* signal the waiting thread */
						pthread_cond_signal(&syncCond);
					}
					pthread_mutex_unlock(&condMutex);
				}
			}

			/* allows another thread to write to this backbuffer */
			void writeBackBuffer(IN_FRAME* value, int bbSeqNo) {
				bbMutex.enter();
				*(this->backBuffer) = *(value);
				this->bbSeqNo = bbSeqNo;
				bbMutex.leave();
			}

			void swapBuffers() {
				IN_FRAME* aux;

				if (POSITION_IN_PIPE == FIRST) {
					/* first stage doesn't need double buffered input */
					return;
				}
				bbMutex.enter();
				aux = frontBuffer;
				frontBuffer = backBuffer;
				backBuffer = aux;
				crtSeqNo = bbSeqNo;
				bbMutex.leave();
			}

			void waitNewInput() {
				if (POSITION_IN_PIPE == FIRST) {
					/*first stage never waits for previous stage */
					return;
				}

				pthread_mutex_lock(&condMutex);
				if (!freshInput) {
					wait4Input = true;
				} else {
					freshInput = false;
				}
				while (wait4Input && !stopFlag) {
					pthread_cond_wait(&syncCond, &condMutex);
				}
				pthread_mutex_unlock(&condMutex);
			}

			/*write the content of the out frame to the back buffer of the next stage */
			void write2Next() {
				/* trick to implement compile time polymorphism */
				write2Next(DummyType<POSITION_IN_PIPE != LAST> ());
			}

			/* deletes all the nodes in the queue placed in front of that which has the same
			 * seq number as the front Buffer;
			 */
			void flushQueue() {
				if (!(SameType<Q_TYPE, NullList>::value::value)) {
					while (readBpQueue->front().seqNo != this->crtSeqNo) {
						readBpQueue->pop(); /* implicitly calls the destructor for the removed element */
					}
				}
			}
	};

	/**
	 * @short the list of types provided to the pipe template
         *
	* NEXT should fit TypeList<br>
	* ELEM is the content
	 */
	template<typename ELEM, typename NEXT = NullList>
	class TypeList {
		public:
			typedef ELEM val;
			typedef NEXT next;
			enum {
				length = 1 + NEXT::length
			};
	};

	/**
	 * @short number of elements in the type list
	 */
	template<typename LIST>
	struct length {
			enum {
				value = LIST::length
			};
	};

	/**
	 * @short class used to iterate through the type list
	 */
	template<typename LIST, int n>
	class NthElem {
		public:
typedef			typename NthElem<typename LIST::next, n - 1>::type_name type_name;
		};

		/**
		 * @short class used to iterate through the type list
		 */
		template<typename LIST>
		class NthElem<LIST, 0> {
			public:
			typedef typename LIST::val type_name;
		};

		/**
		 * @short retrieve the q_type as the out_frame_type of the stage at the given offset
		 */
		template<typename LIST, int step, int q_from>
		class GetQType {
			public:
			typedef typename NthElem<LIST, step - q_from>::type_name::out_frame_type q_type;
		};

		/**
		 * @short retrieve the q_type as the out_frame_type of the stage at the given offset
		 */
		template<typename LIST, int step>
		class GetQType<LIST, step, 0> {
			public:
			typedef NullList q_type;
		};

		/**
		 * @short classes used to build instances of the stages in the pipeline and the links between them:
		 * interface for a node in the pipeline
		 */
		class NodeBase {
			public:
			void* nextNode; /* next node in the PipeNode list (the list used when building the pipe) */
			int qFromDelta; /* offset to the stage of origin for the bypass q */
			int pipe_step; /* index in the list used when building the pipe */

			virtual void setWriteQ(void* writeQ)=0;
			virtual void* getReadQ()=0;
			virtual void start()=0;
			virtual void stop()=0;
			virtual void join()=0;
		};

		/**
		 * @short classes used to build instances of the stages in the pipeline and the links between them
		 */
		template<typename LIST, int n, int step, typename NEXT_STAGE>
		class PipeNodes: public NodeBase {
			public:
			typedef typename LIST::val stage_types;
			typedef GenericStage<typename LIST::next::val::out_frame_type, typename stage_types::out_frame_type, NEXT_STAGE,
			typename GetQType<LIST, step, stage_types::q_from>::q_type, typename stage_types::computation> crt_stage_type;

			crt_stage_type value;

			PipeNodes(NEXT_STAGE* next) {
				value.setNext(next);
				qFromDelta = step - stage_types::q_from;
				pipe_step = step;
				nextNode = new PipeNodes<typename LIST::next, n + 1, step - 1, crt_stage_type> (&value);
			}

			void* getReadQ() {
				return value.getReadQ();
			}

			void setWriteQ(void* writeQ) {
				value.setWriteQ((std::queue<Queueable<typename stage_types::out_frame_type> >*)writeQ);
			}

			void start() {
				value.start();
			}

			void stop() {
				value.stop();
			}

			void join() {
				value.join();
			}
		};

		/**
		 * @short classes used to build instances of the stages in the pipeline and the links between them:
		 * first stage in the pipeline
		 */
		template<typename LIST, int n, typename NEXT_STAGE>
		class PipeNodes<LIST, n, 1, NEXT_STAGE> : public NodeBase {
			public:
			typedef typename LIST::val stage_types;
			/* first stage thus no bypass Q to read from */
			typedef GenericStage<NullList, typename stage_types::out_frame_type, NEXT_STAGE,
			NullList, typename stage_types::computation, FIRST> crt_stage_type;

			crt_stage_type value;

			PipeNodes(NEXT_STAGE* next) {
				value.setNext(next);
				nextNode = NULL;
				qFromDelta = 1;
				pipe_step = 1;
			}

			void* getReadQ() {
				return value.getReadQ();
			}

			void setWriteQ(void* writeQ) {
				value.setWriteQ((std::queue<Queueable<typename stage_types::out_frame_type> >*)writeQ);
			}

			void start() {
				value.start();
			}

			void stop() {
				value.stop();
			}

			void join() {
				value.join();
			}
		};

		/**
		 * @short classes used to build instances of the stages in the pipeline and the links between them:
		 * last stage in the pipeline
		 */
		template<typename LIST, int step>
		class PipeNodes<LIST, 0, step, NullList> : public NodeBase {
			public:
			typedef typename LIST::val stage_types;
			typedef GenericStage<typename LIST::next::val::out_frame_type, typename stage_types::out_frame_type, NullList,
			typename GetQType<LIST, step, stage_types::q_from>::q_type, typename stage_types::computation, LAST> crt_stage_type;

			crt_stage_type value;

			PipeNodes() {
				qFromDelta = step - stage_types::q_from;
				pipe_step = step;
				nextNode = new PipeNodes<typename LIST::next, 1, step - 1, crt_stage_type> (&value);
			}

			void* getReadQ() {
				return value.getReadQ();
			}

			void setWriteQ(void* writeQ) {
				value.setWriteQ((std::queue<Queueable<typename stage_types::out_frame_type> >*)writeQ);
			}

			void start() {
				value.start();
			}

			void stop() {
				value.stop();
			}

			void join() {
				value.join();
			}
		};

		/**
		 * @short builds the nodes in the pipeline and the links between the stages
		 */
		template<typename LIST>
		class Pipe {
			private:
			PipeNodes<LIST, 0, length<LIST>::value, NullList>
			pipeNodes;

			public:
			/**
			 * @short class constructor
			 */
			Pipe() {
				/* set the sources for the bypass queues */
				NodeBase* crt = &pipeNodes;
				while (crt->nextNode != NULL) {
					if (crt->qFromDelta != crt->pipe_step) {
						NodeBase* ptr = crt;
						do {
							ptr = (NodeBase*) (ptr->nextNode);
							crt->qFromDelta--;
						}while (crt->qFromDelta> 0);
						ptr ->setWriteQ(crt->getReadQ());
					}
					crt = (NodeBase*) (crt->nextNode);
				}
			}

			/** @short start all the stages in the pipeline */
			void start() {
				NodeBase* crt = &pipeNodes;
				while (crt->nextNode != NULL) {
					crt->start();
					crt = (NodeBase*) (crt->nextNode);
				}
				crt->start();
			}

			/** @short stop all the stages in the pipeline */
			void stop() {
				NodeBase* crt = &pipeNodes;
				while (crt->nextNode != NULL) {
					crt->stop();
					crt = (NodeBase*) (crt->nextNode);
				}
				crt->stop();
			}

			/** @short wait for the threads to stop */
			void join() {
				NodeBase* crt = &pipeNodes;
				while (crt->nextNode != NULL) {
					crt->join();
					crt = (NodeBase*) (crt->nextNode);
				}
				crt->join();
			}

		};
	}
/**
* @example pipeliner_example_00.cc */

#endif /* PIPELINE_H_ */
