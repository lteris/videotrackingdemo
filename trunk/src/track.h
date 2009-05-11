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
	extern std::string param_buffer_in_hostname;
	extern std::string param_buffer_in_resource;
	extern std::string param_buffer_in_entity;
	extern int param_buffer_in_port;
	extern bool param_buffer_in_jpg;
	extern std::string param_buffer_out_hostname;
	extern std::string param_buffer_out_resource;
	extern std::string param_buffer_out_entity;
	extern int param_buffer_out_port;
	extern bool param_buffer_out_jpg;
	extern int param_out_jpeg_quality;

	extern std::string param_buffer_detection_hostname;
	extern std::string param_buffer_detection_resource;
	extern std::string param_buffer_detection_entity;
	extern int param_buffer_detection_port;
	extern bool param_buffer_detection_jpg;
	extern int param_detection_jpeg_quality;
	extern bool param_display_detection;

	extern double param_gngt_target;
	extern double param_gngt_first_learning_rate;
	extern double param_gngt_second_learning_rate;
	extern int param_gngt_edge_age_max;
	extern double param_gngt_variance_max;
	extern double param_gngt_length_max;
	extern int param_nb_epochs_per_frame;
	extern int param_pen_thickness;

	extern bool param_morph;
	extern int param_morph_radius;

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
			/* post image buffers */
			bkbd::Image dest;
			jpeg::Compress compressor;
			bkbd::JPEG outJpegBuf;
			/* get image buffers */
			bkbd::JPEG inJpegBuf;
			jpeg::Decompress decompressor;
			struct timeval timestamp;
		public:
			ServerConnection();
			virtual ~ServerConnection();
			/* retrieve an image from the input entity */
			void getImage(bkbd::Image* frame);
			/* post outFrame to the output entity */
			void postImg(ImageRGB24& outFrame);

			/* test method */
			void postImg(ImageBool& outFrame);
	};

	extern ServerConnection* serverConn;

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
	class ImageFeeder_RGB24 {
		private:
			bkbd::Image buffer;
		public:
			void operator()(void* dummy, ImageRGB24*& outFrame);
	};

	/*-----------------------------------------------------------------------------------------*/
	/* posts images to the server                                                              */
	/*-----------------------------------------------------------------------------------------*/

	class ImagePoster {
		public:
			void operator()(ImageRGB24* image, void* dummy);
			void operator()(ImageBool* image, void* dummy);
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
			ImageGRAY8 buffer;
			mirage::img::ForegroundDetection<DetectionParams, 0> algo;
			bool firstImg;

			void convert2Gray(ImageRGB24* inFrame);
			void getForeground(ImageBool*& outFrame);
		public:
			ForeGround() {
				firstImg = true;
			}
			void operator()(ImageRGB24* inFrame, ImageBool*& outFrame){
				convert2Gray(inFrame);
				getForeground(outFrame);
			}
	};

	/*-----------------------------------------------------------------------------------------*/
	/* various transformations on the image                                                    */
	/*-----------------------------------------------------------------------------------------*/
	class Morpho_Contour {
		private:
			typedef mirage::img::Coding<mirage::morph::MaskValue>::Frame
					Element;
			Element element;
			ImageBool buffer;

			void morphoMath(ImageBool* inFrame);
			void getContour(ImageBool& outFrame);
		public:
			Morpho_Contour() {
				mirage::morph::Mask::Disk(element, param_morph_radius);
			}
			void operator()(ImageBool* inFrame, ImageBool*& outFrame){
				morphoMath(inFrame);
				getContour(*outFrame);
			}
	};

	/*-----------------------------------------------------------------------------------------*/
	/* applies gngt on the contour retrieved by the contour detection class                    */
	/*-----------------------------------------------------------------------------------------*/
	class GNGT_Draw {
		private:
			GNG_T algo;
			std::vector<Input> examples;
			std::map<int, mirage::colorspace::RGB_24> colors;
			mirage::img::Coordinate pen, half_pen;
			LABELIZER labelizer;

			mirage::colorspace::RGB_24 getRandomColor();

			void feedGNGT(ImageBool& contour);
			void labelize();
			void draw(ImageRGB24*& outFrame, ImageRGB24& fromQueue);
		public:
			GNGT_Draw();
			void operator()(ImageBool* contour, ImageRGB24*& outFrame,
					ImageRGB24& fromQueue) {
				feedGNGT(*contour);
				labelize();
				draw(outFrame, fromQueue);
			}

	};
}

#endif /* TRACK_H_ */
