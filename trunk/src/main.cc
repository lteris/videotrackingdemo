#include "track.h"
#include "pipeline.h"

using namespace pipeline;

int main() {
	typedef Types4Stage<track::ImageRGB24, track::ImagePoster, 0> stage_poster;
	typedef Types4Stage<track::ImageRGB24, track::GNGT_Draw, 1> stage_gngt_draw;
	typedef Types4Stage<track::ImageBool, track::Morpho_Contour, 0> stage_morpho_contour;
	typedef Types4Stage<track::ImageBool, track::ForeGround, 0>	stage_foreground;
	typedef Types4Stage<track::ImageRGB24, track::ImageFeeder_RGB24, 0> stage_feeder;

	typedef TypeList<stage_feeder> list_1;
	typedef TypeList<stage_foreground, list_1> list_2;
	typedef TypeList<stage_morpho_contour, list_2> list_3;
	typedef TypeList<stage_gngt_draw, list_3> list_4;
	typedef TypeList<stage_poster, list_4> type_list;

	track::ParameterParser params("toto.demo");
	Pipe<type_list> pipe;

	pipe.start();

	//	ost::Thread::sleep(1000000);
	while (true)
		;

	pipe.stop();
	pipe.join();
}

