#include "track.h"
#include "pipeline.h"
#include "test.h"

using namespace pipeline;

int main() {
	typedef Types4Stage<ImageRGB24, ImagePoster, 0> stage_poster;
	typedef Types4Stage<ImageRGB24, Draw, 2> stage_draw;
	typedef Types4Stage<GNGTGraph, GNGT, 0> stage_gngt;
	typedef Types4Stage<ImageBool, Contour, 0> stage_contour;
	typedef Types4Stage<ImageBool, MorphoMath, 0> stage_morpho;
	typedef Types4Stage<ImageBool, ForeGround, 0> stage_foreground;
	typedef Types4Stage<ImageGRAY8, Convert2Gray, 0> stage_get_gray;
	typedef Types4Stage<ImageRGB24, Convert2RGB, 0> stage_get_rgb;
	typedef Types4Stage<bkbd::Image, ImageFeeder, 0> stage_feeder;

	typedef TypeList<stage_feeder> list_1;
	typedef TypeList<stage_get_rgb, list_1> list_2;
	typedef TypeList<stage_get_gray, list_2> list_3;
	typedef TypeList<stage_foreground, list_3> list_4;
	typedef TypeList<stage_morpho, list_4> list_5;
	typedef TypeList<stage_contour, list_5> list_6;
	typedef TypeList<stage_gngt, list_6> list_7;
	typedef TypeList<stage_draw, list_7> list_8;
	typedef TypeList<stage_poster, list_8> type_list;

	Pipe<type_list, ParameterParser> pipe("config.xml");

	pipe.start();

	ost::Thread::sleep(2000);

	pipe.stop();
	pipe.join();
}

