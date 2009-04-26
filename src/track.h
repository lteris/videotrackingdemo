#ifndef TRACK_H_
#define TRACK_H_

#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <string>
#include <mirage.h>
#include <bkbd.h>
#include <cc++/thread.h>
#include <vq.h>
#include <gngtlib.h>
#include <vector>
#include <algorithm>
#include <map>
#include <string>
#include <queue>
#include "Foreground.h"

namespace track {

	int param_inter_frame_delay;
	std::string param_buffer_in_hostname;
	std::string param_buffer_in_resource;
	std::string param_buffer_in_entity;
	int param_buffer_in_port;
	bool param_buffer_in_jpg;
	std::string param_buffer_out_hostname;
	std::string param_buffer_out_resource;
	std::string param_buffer_out_entity;
	int param_buffer_out_port;
	bool param_buffer_out_jpg;
	int param_out_jpeg_quality;
	std::string param_buffer_detection_hostname;
	std::string param_buffer_detection_resource;
	std::string param_buffer_detection_entity;
	int param_buffer_detection_port;
	bool param_buffer_detection_jpg;
	int param_detection_jpeg_quality;
	bool param_display_detection;
	double param_color_min_norm2;
	double param_color_min_norm2_margin;
	double param_color_min_cosine;
	int param_nb_background_samples;
	double param_background_coef;
	double param_gngt_target;
	double param_gngt_first_learning_rate;
	double param_gngt_second_learning_rate;
	int param_gngt_edge_age_max;
	double param_gngt_variance_max;
	double param_gngt_length_max;
	int param_nb_epochs_per_frame;
	int param_pen_thickness;
	int param_left_margin;
	int param_right_margin;
	int param_top_margin;
	int param_bottom_margin;
	bool param_morph;
	int param_morph_radius;

	typedef mirage::img::Coding<mirage::colorspace::GRAY_8>::Frame ImageGRAY8;
	typedef mirage::img::Coding<mirage::colorspace::RGB_24>::Frame ImageRGB24;
	typedef mirage::img::Coding<mirage::colorspace::RGB<double> >::Frame
			ImageRGBf;
	typedef mirage::img::Coding<bool>::Frame ImageBool;
	typedef mirage::img::Coding<mirage::morph::MaskValue>::Frame Element;
	typedef vq::Vector<2, double> Input;

	class ParamGNG_T: public gngt::ParamVQ<Input> {
		public:
			static int AgeMax(void) {
				return param_gngt_edge_age_max;
			}
			static double FirstLearningRate(void) {
				return param_gngt_first_learning_rate;
			}
			static double SecondLearningRate(void) {
				return param_gngt_second_learning_rate;
			}
			static double Target(void) {
				return param_gngt_target;
			}
			static double Lambda(void) {
				return .9;
			}
	};

	typedef gngt::Network<ParamGNG_T> GNG_T;
	typedef vq::Labelizer<GNG_T> LABELIZER;

	/*-----------------------------------------------------------------------------------------*/
	/* retrieve images from bkbd server _ manages server connection                            */
	/* IN_FRAME : type of images retrieved from the server									   */
	/* OUT_FRAME: tye of output produced by last stage i.e. output posted on the server        */
	/*-----------------------------------------------------------------------------------------*/
	class ServerConnection {
		private:
			bkbd::Repository<bkbd::JPEG>::Client* jpgClientIn;
			bkbd::Repository<bkbd::JPEG>::Client* jpgClientOut;
			bkbd::Repository<bkbd::JPEG>::Client* jpgClientDetection;
			bkbd::Repository<bkbd::Image>::Client* rawClientIn;
			bkbd::Repository<bkbd::Image>::Client* rawClientOut;
			bkbd::Repository<bkbd::Image>::Client* rawClientDetection;
		public:
			ServerConnection();
			virtual ~ServerConnection();
			/* retrieve an image from the input entity */
			void getImage(bkbd::Image* frame);
			/* post outFrame to the output entity */
			void postImg(ImageRGB24& outFrame);
	};

	ServerConnection* serverConn;

	/*-----------------------------------------------------------------------------------------*/
	/* loads parameters form file or default parameters                                        */
	/*-----------------------------------------------------------------------------------------*/
	class ParameterParser {
		private:
#define Mode(jpg)((jpg)?"jpeg":"rgb24")
			void loadDefaultParameters();
		public:
			/* default parameters */
			ParameterParser();
			/* load parameters from file */
			ParameterParser(const std::string& confFile);
			void saveParameters(const std::string& filename);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* retrieves images from the server and feeds them to the next stage                       */
	/*-----------------------------------------------------------------------------------------*/
	class ImageFeeder {
		public:
			void operator()(void* dummy, bkbd::Image*& outFrame);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* posts images to the server                                                              */
	/*-----------------------------------------------------------------------------------------*/

	class ImagePoster {
		public:
			void operator()(ImageRGB24* image, void* dummy);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* obtain mirage rgb24 image image from bkbd::Image                                          */
	/*-----------------------------------------------------------------------------------------*/
	class Convert2RGB {
		public:
			void operator()(bkbd::Image* inFrame, ImageRGB24*& outFrame);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* obtain grayscale image from rgb24                                                       */
	/*-----------------------------------------------------------------------------------------*/
	class Convert2Gray {
		public:
			void operator()(ImageRGB24* inFrame, ImageGRAY8*& outFrame);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* identify the background pixels                                                          */
	/*-----------------------------------------------------------------------------------------*/
	class ForeGround {
		private:
			class DetectionParams {
				public:
					inline static int stdMin(void) {
						return 3;
					}
					inline static int stdMax(void) {
						return 20;
					}
					inline static double k(void) {
						return 2.5;
					}
					inline static double shadowRatio(void) {
						return 0.55;
					}
					inline static double shadowVar(void) {
						return 0.075;
					}
					inline static int shadowWindow(void) {
						return 1;
					}
					inline static int morphEro(void) {
						return 1;
					}
					inline static int morphDil(void) {
						return 2;
					}
			};
			mirage::img::ForegroundDetection<DetectionParams, 0> algo;
			bool firstImg;
		public:
			ForeGround() {
				firstImg = true;
			}
			void operator()(ImageGRAY8* inFrame, ImageBool*& outFrame);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* various transformations on the image                                                    */
	/*-----------------------------------------------------------------------------------------*/
	class MorphoMath {
		private:
			typedef mirage::img::Coding<mirage::morph::MaskValue>::Frame
					Element;
			Element element;
		public:
			void operator()(ImageBool* inFrame, ImageBool*& outFrame);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* retrieves contour to be fed to the GNGT module							               */
	/*-----------------------------------------------------------------------------------------*/
	class Contour {
		public:
			void operator()(ImageBool* inFrame, ImageBool*& outFrame);

	};

	/*-----------------------------------------------------------------------------------------*/
	/* applies gngt on the contour retrieved by the contour detection class                    */
	/*-----------------------------------------------------------------------------------------*/
	class GNGT {
		private:
			GNG_T algo;
			std::vector<Input> examples;

			void feedGNGT(ImageBool& contour);
			void labelize(LABELIZER& labelizer);
		public:
			void operator()(ImageBool* contour, LABELIZER*& labelizer) {
				feedGNGT(*contour);
				labelize(*labelizer);
			}

	};

	/*-----------------------------------------------------------------------------------------*/
	/* draws the graph over the image received from the conversion stage                       */
	/*-----------------------------------------------------------------------------------------*/
	class Draw {
		private:
			std::map<int, mirage::colorspace::RGB_24> colors;
			ImageRGB24::point_type pen, half_pen;

			mirage::colorspace::RGB_24 getRandomColor();
		public:
			Draw();
			void operator()(LABELIZER* inFrame, ImageRGB24*& outFrame,
					ImageRGB24& fromQueue);
	};

}

#endif /* TRACK_H_ */