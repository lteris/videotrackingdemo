/*
 * gngt_track.h
 *
 *  Created on: May 12, 2009
 *      Author: liviu
 */

#ifndef GNGT_TRACK_H_
#define GNGT_TRACK_H_

#include "track.h"
#include <pipeline.h>

namespace gngt_track {

	/**
	 * builds the pipeline
	 */
	template<typename SERVER>
	class VideoTracker {
		private:
			typedef pipeline::Types4Stage<typename track<SERVER>::ImageRGB24,
					typename track<SERVER>::ImagePoster, 0> stage_poster;
			typedef pipeline::Types4Stage<typename track<SERVER>::ImageRGB24,
					typename track<SERVER>::GNGT_Draw, 1> stage_gngt_draw;
			typedef pipeline::Types4Stage<typename track<SERVER>::ImageBool,
					typename track<SERVER>::Morpho_Contour, 0>
					stage_morpho_contour;
			typedef pipeline::Types4Stage<typename track<SERVER>::ImageBool,
					typename track<SERVER>::ForeGround, 0> stage_foreground;
			typedef pipeline::Types4Stage<typename track<SERVER>::ImageRGB24,
					typename track<SERVER>::ImageFeeder_RGB24, 0> stage_feeder;

			typedef pipeline::TypeList<stage_feeder> list_1;
			typedef pipeline::TypeList<stage_foreground, list_1> list_2;
			typedef pipeline::TypeList<stage_morpho_contour, list_2> list_3;
			typedef pipeline::TypeList<stage_gngt_draw, list_3> list_4;
			typedef pipeline::TypeList<stage_poster, list_4> type_list;

			pipeline::Pipe<type_list>* pipe;
			bool srvAlloced;
		public:
			/**
			 * builds the pipeline
			 * @param gngtConf configuration file containing the gngt parameters
			 * @param srvConf server configuration file
			 */
			VideoTracker(const std::string& gngtConf,
					const std::string& srvConf) {
				track<SERVER>::parameterParser(gngtConf);
				track<SERVER>::serverConn() = new SERVER(srvConf);
				srvAlloced = true;
				pipe = new pipeline::Pipe<type_list>();
			}

			/**
			 * builds the pipeline
			 * @param gngtConf configuration file containing the gngt parameters
			 * @param srv externally built server
			 */
			VideoTracker(const std::string& gngtConf, SERVER* srv) {
				track<SERVER>::serverConn() = srv;
				track<SERVER>::parameterParser(gngtConf);
				pipe = new pipeline::Pipe<type_list>();
				srvAlloced = false;
			}

			~VideoTracker() {
				delete pipe;
				/* delete de server unless allocated externally */
				if (srvAlloced) {
					delete track<SERVER>::serverConn();
				}
			}

			/**
			 * commence pipeline operation - all stages are blocked until they receive some input
			 */
			void start() {
				pipe->start();
			}

			/**
			 * stop all stages in the pipeline
			 */
			void stop() {
				pipe->stop();
			}

			/**
			 * wait for the stages to stop
			 */
			void join() {
				pipe->join();
			}
	};

	template <typename SERVER>
	class BackgroundTracker {
			private:
				typedef pipeline::Types4Stage<typename track<SERVER>::ImageRGB24,
						typename track<SERVER>::ImagePoster, 0> stage_poster;
				typedef pipeline::Types4Stage<typename track<SERVER>::ImageBool,
						typename track<SERVER>::ForeGround, 0> stage_foreground;
				typedef pipeline::Types4Stage<typename track<SERVER>::ImageRGB24,
						typename track<SERVER>::ImageFeeder_RGB24, 0> stage_feeder;

				typedef pipeline::TypeList<stage_feeder> list_1;
				typedef pipeline::TypeList<stage_foreground, list_1> list_2;
				typedef pipeline::TypeList<stage_poster, list_2> type_list;

				pipeline::Pipe<type_list>* pipe;
				bool srvAlloced;
			public:
				/**
				 * builds the pipeline
				 * @param gngtConf configuration file containing the gngt parameters
				 * @param srvConf server configuration file
				 */
				BackgroundTracker(const std::string& gngtConf,
						const std::string& srvConf) {
					track<SERVER>::parameterParser(gngtConf);
					track<SERVER>::serverConn() = new SERVER(srvConf);
					srvAlloced = true;
					pipe = new pipeline::Pipe<type_list>();
				}

				/**
				 * builds the pipeline
				 * @param gngtConf configuration file containing the gngt parameters
				 * @param srv externally built server
				 */
				BackgroundTracker(const std::string& gngtConf, SERVER* srv) {
					track<SERVER>::serverConn() = srv;
					track<SERVER>::parameterParser(gngtConf);
					pipe = new pipeline::Pipe<type_list>();
					srvAlloced = false;
				}

				~BackgroundTracker() {
					delete pipe;
					/* delete de server unless allocated externally */
					if (srvAlloced) {
						delete track<SERVER>::serverConn();
					}
				}

				/**
				 * commence pipeline operation - all stages are blocked until they receive some input
				 */
				void start() {
					pipe->start();
				}

				/**
				 * stop all stages in the pipeline
				 */
				void stop() {
					pipe->stop();
				}

				/**
				 * wait for the stages to stop
				 */
				void join() {
					pipe->join();
				}
		};


}

#endif /* GNGT_TRACK_H_ */
