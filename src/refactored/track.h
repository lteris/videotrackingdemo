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
#include <foreground.h>

#define KWD_TARGET            "GNGT_TARGET"
#define KWD_AGE               "GNGT_EDGE_AGE"
#define KWD_LEARNING_RATES    "GNGT_LEARNING_RATES"
#define KWD_VARIANCE          "GNGT_MAX_VARIANCE"
#define KWD_LENGTH            "GNGT_MAX_LENGTH"
#define KWD_EPOCHS            "GNGT_EPOCHS_PER_FRAME"
#define KWD_PEN_THICKNESS     "PEN_THICKNESS"
#define KWD_MORPH             "MORPHOMATH"

template<typename SERVER>
class track {
		static double& param_gngt_target() {
			static double param_gngt_target;
			return param_gngt_target;
		}
		static double& param_gngt_first_learning_rate() {
			static double param_gngt_first_learning_rate;
			return param_gngt_first_learning_rate;
		}
		static double& param_gngt_second_learning_rate() {
			static double param_gngt_second_learning_rate;
			return param_gngt_second_learning_rate;
		}

		static int& param_gngt_edge_age_max() {
			static int param_gngt_edge_age_max;
			return param_gngt_edge_age_max;
		}

		static double& param_gngt_variance_max() {
			static double param_gngt_variance_max;
			return param_gngt_variance_max;
		}
		static double& param_gngt_length_max() {
			static double param_gngt_length_max;
			return param_gngt_length_max;
		}

		static int& param_nb_epochs_per_frame() {
			static int param_nb_epochs_per_frame;
			return param_nb_epochs_per_frame;
		}

		static int& param_pen_thickness() {
			static int param_pen_thickness;
			return param_pen_thickness;
		}

		static bool& param_morph() {
			static bool param_morph;
			return param_morph;
		}

		static int& param_morph_radius() {
			static int param_morph_radius;
			return param_morph_radius;
		}

	public:
		typedef mirage::img::Coding<mirage::colorspace::GRAY_8>::Frame
				ImageGRAY8;
		typedef mirage::img::Coding<mirage::colorspace::RGB_24>::Frame
				ImageRGB24;
		typedef mirage::img::Coding<mirage::colorspace::RGB<double> >::Frame
				ImageRGBf;
		typedef mirage::img::Coding<bool>::Frame ImageBool;
		typedef mirage::img::Coding<mirage::morph::MaskValue>::Frame Element;
		typedef vq::Vector<2, double> Input;
	public:
		class ParamGNG_T: public gngt::ParamVQ<Input> {
			public:
				static int AgeMax(void) {
					return param_gngt_edge_age_max();
				}
				static double FirstLearningRate(void) {
					return param_gngt_first_learning_rate();
				}
				static double SecondLearningRate(void) {
					return param_gngt_second_learning_rate();
				}
				static double Target(void) {
					return param_gngt_target();
				}
				static double Lambda(void) {
					return .9;
				}
		};

		typedef gngt::Network<ParamGNG_T> GNG_T;
		typedef vq::Labelizer<GNG_T> LABELIZER;

		static SERVER*& serverConn() {
			static SERVER* serverConn;
			return serverConn;
		}

#define Mode(jpg)((jpg)?"jpeg":"rgb24")
		static void loadDefaultParameters() {
			param_gngt_target() = 30;
			param_gngt_first_learning_rate() = .05;
			param_gngt_second_learning_rate() = .005;
			param_gngt_edge_age_max() = 20;
			param_gngt_variance_max() = 50;
			param_gngt_length_max() = 20;
			param_nb_epochs_per_frame() = 10;
			param_pen_thickness() = 3;
			param_morph() = false;
			param_morph_radius() = 3;
		}

		/**
		 * set default parameters
		 */
		static void parameterParser() {
			loadDefaultParameters();
		}

		/**
		 * get the parameters from the input file
		 */
		static void parameterParser(const std::string& filename) {
			std::ifstream file;
			std::string comment;
			std::string kwd;
			std::string mode;
			char first;

			loadDefaultParameters();

			file.open(filename.c_str());
			if (!file) {
				std::cerr << "Cannot open file \"" << filename
						<< "\" for loading configuration." << std::endl;
				::exit(1);
			}

			file >> std::ws;
			file.get(first);
			while (!file.eof()) {

				file.putback(first);
				if (first == '#')
					std::getline(file, comment, '\n');
				else {
					file >> kwd;
					if (kwd == KWD_TARGET) {
						file >> param_gngt_target();
					} else if (kwd == KWD_AGE)
						file >> param_gngt_edge_age_max();
					else if (kwd == KWD_LEARNING_RATES)
						file >> param_gngt_first_learning_rate()
								>> param_gngt_second_learning_rate();
					else if (kwd == KWD_VARIANCE)
						file >> param_gngt_variance_max();
					else if (kwd == KWD_LENGTH)
						file >> param_gngt_length_max();
					else if (kwd == KWD_EPOCHS)
						file >> param_nb_epochs_per_frame();
					else if (kwd == KWD_PEN_THICKNESS) {
						file >> param_pen_thickness();
					} else if (kwd == KWD_MORPH) {
						param_morph() = true;
						file >> param_morph_radius();
					} else {
						std::cerr << "Warning : Keyword \"" << kwd
								<< "\" is ignored, line skipped." << std::endl;
						std::getline(file, comment, '\n');
					}
				}

				file >> std::ws;
				file.get(first);
			}

			file.close();
		}

		/**
		 * save current parameters to file
		 */
		static void saveParameters(const std::string& filename) {
			std::ofstream file;

			file.open(filename.c_str());
			if (!file) {
				std::cerr << "Cannot open file \"" << filename
						<< "\" for saving configuration." << std::endl;
				::exit(1);
			}

			file << "# GNG-T Parameters : target, maximum edge age,"
					<< std::endl << "# winner and second learning rates."
					<< std::endl << KWD_TARGET << ' ' << param_gngt_target()
					<< std::endl << KWD_AGE << ' ' << param_gngt_edge_age_max()
					<< std::endl << KWD_LEARNING_RATES << ' '
					<< param_gngt_first_learning_rate() << ' '
					<< param_gngt_second_learning_rate() << std::endl
					<< "# Number of GNG-T epochs per frame (a supplementary one will be actually used)."
					<< std::endl << KWD_EPOCHS << ' '
					<< param_nb_epochs_per_frame() << std::endl
					<< "# Maximal variance for a non noisy node." << std::endl
					<< KWD_VARIANCE << ' ' << param_gngt_variance_max()
					<< std::endl << "# Maximal edge length inside a polygon."
					<< std::endl << KWD_LENGTH << ' '
					<< param_gngt_length_max() << std::endl << std::endl
					<< "# Drawing paremeters" << std::endl << KWD_PEN_THICKNESS
					<< ' ' << param_pen_thickness() << std::endl << std::endl
					<< "# Uncomment the following line to enable morpho-matematical"
					<< std::endl
					<< "# detection cleaning. Number is the mask radius."
					<< std::endl;
			if (!param_morph())
				file << "# ";
			file << KWD_MORPH << ' ' << param_morph_radius() << std::endl
					<< std::endl << std::endl;

			file.close();
		}

		/**
		 * retrieves bkbd::Image from the server, convert them to ImageRGB24
		 */
		class ImageFeeder_RGB24 {
			private:
				bkbd::Image buffer;
			public:
				/**
				 * @param dummy the first stage does not have an input value - value used so that the class
				 * follows the COMPUTE concept in pipeliner
				 * @param outFrame - output produced by this stage
				 */
				void operator()(void* dummy, ImageRGB24*& outFrame) {
					serverConn()->getImage(&buffer);
					outFrame->resize(mirage::img::Coordinate(buffer.width,
							buffer.height),
							(mirage::colorspace::RGB_24*) (buffer.data));
				}
		};

		/**
		 * posts ImageRGB24 to the server
		 */
		class ImagePoster {
			public:
				/**
				 * @param image the image to post
				 * @param dummy the last stage does not have an out put value - value used so that the class
				 * follows the COMPUTE concept in pipeliner
				 */
				void operator()(ImageRGB24* image, void* dummy) {
					serverConn()->postImg(*image);
				}

				void operator()(ImageBool* image, void* dummy) {
					serverConn()->postImg(*image);
				}
		};

		/**
		 * identify the foreground pixels
		 */
		class ForeGround {
			private:
				/**
				 * parameters for the foreground detection algorithm
				 */
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

				void convert2Gray(ImageRGB24* inFrame) {
					buffer.resize(inFrame->_dimension);
					mirage::algo::UnaryOp<ImageRGB24, ImageGRAY8,
							mirage::colorspace::RGBToGray<
									ImageRGB24::value_type,
									ImageGray8::value_type> >(*inFrame, buffer);
				}
				void getForeground(ImageBool*& outFrame) {
					if (firstImg) {
						algo.initForegroundDetection(buffer);
						firstImg = false;
					}
					algo.SigmaDeltaModified(buffer);
					algo.getMask(*outFrame);
				}
			public:
				ForeGround() {
					firstImg = true;
				}

				/**
				 * identify the foreground pixels
				 * @param inFrame - input image
				 * @param outFrame - image containing the foreground pixels
				 */
				void operator()(ImageRGB24* inFrame, ImageBool*& outFrame) {
					convert2Gray(inFrame);
					getForeground(outFrame);
				}
		};

		/*
		 * apply various transformations on the image and get the contour of the foreground
		 */
		class Morpho_Contour {
			private:
				typedef mirage::img::Coding<mirage::morph::MaskValue>::Frame
						Element;
				Element element;
				ImageBool buffer;

				void morphoMath(ImageBool* inFrame) {
					ImageBool::pixel_type pix1, pix2, pix_end;
					ImageBool::point_type pos, offset;
					bool n, s, e, w;

					buffer.resize(inFrame->_dimension);
					if (!param_morph()) {
						for (pix1 = inFrame->begin(), pix_end = inFrame->end(), pix2
								= buffer.begin(); pix1 != pix_end; ++pix1, ++pix2) {
							if (*pix1) {
								pos = !pix2 + offset(1, 0);
								w = !((*inFrame)(pos));
								pos = !pix2 + offset(-1, 0);
								e = !((*inFrame)(pos));
								pos = !pix2 + offset(0, 1);
								s = !((*inFrame)(pos));
								pos = !pix2 + offset(0, -1);
								n = !((*inFrame)(pos));

								*pix2 = !((n && s) || (e && w));
							} else {
								*pix2 = false;
							}
						}
					} else {
						mirage::morph::Format<ImageBool, ImageBool, 0>::Opening(
								*inFrame, element, buffer);
					}
				}

				void getContour(ImageBool& outFrame) {
					ImageBool::pixel_type pix1, pix2, pix_end;
					ImageBool::point_type pos, dimension, offset;

					outFrame.resize(buffer._dimension);
					for (pix1 = buffer.begin(), pix_end = buffer.end(), pix2
							= outFrame.begin(); pix1 != pix_end; ++pix1, ++pix2) {
						if (*pix1) {
							pos = !pix2 + offset(1, 0);

							if (!(*pix2 = !(buffer(pos)))) {
								pos = !pix2 + offset(0, 1);
								if (!(*pix2 = !(buffer(pos)))) {
									pos = !pix2 + offset(0, -1);
									if (!(*pix2 = !(buffer(pos)))) {
										pos = !pix2 + offset(-1, 0);
										*pix2 = !(buffer(pos));
									}
								}
							}

						} else {
							*pix2 = false;
						}
					}
				}
			public:
				Morpho_Contour() {
					mirage::morph::Mask::Disk(element, param_morph_radius());
				}

				/**
				 * identify the foreground pixels
				 * @param inFrame - input image
				 * @param outFrame - image containing the contour of the foreground
				 */
				void operator()(ImageBool* inFrame, ImageBool*& outFrame) {
					morphoMath(inFrame);
					getContour(*outFrame);
				}
		};

		/*
		 * applies gngt on the contour retrieved by the contour detection class
		 */
		class GNGT_Draw {
			private:
				GNG_T algo;
				std::vector<Input> examples;
				std::map<int, mirage::colorspace::RGB_24> colors;
				mirage::img::Coordinate pen, half_pen;
				LABELIZER labelizer;

				mirage::colorspace::RGB_24 getRandomColor() {
					mirage::colorspace::RGB_24 res;

					res._red = (mirage::colorspace::RGB_24::value_type) (256.0
							* (rand() / (RAND_MAX + 1.0)));
					res._green
							= (mirage::colorspace::RGB_24::value_type) (256.0
									* (rand() / (RAND_MAX + 1.0)));
					res._blue = (mirage::colorspace::RGB_24::value_type) (256.0
							* (rand() / (RAND_MAX + 1.0)));

					return res;
				}

				void feedGNGT(ImageBool& contour) {
					ImageBool::pixel_type pix, pix_end;
					std::vector<Input>::iterator iter, iter_end;
					Input example;
					int epoch;

					/* Getting examples */
					examples.clear();
					for (pix = contour.begin(), pix_end = contour.end(); pix
							!= pix_end; ++pix)
						if (*pix) {
							example[0] = (!pix)[0];
							example[1] = (!pix)[1];
							examples.push_back(example);
						}

					for (epoch = 0; epoch < param_nb_epochs_per_frame(); ++epoch) {
						iter_end = examples.end();
						iter = examples.begin();
						std::random_shuffle(iter, iter_end);

						algo.OpenEpoch(true);
						for (; iter != iter_end; ++iter)
							algo.Submit(*iter);
						algo.CloseEpoch();
					}

					iter_end = examples.end();
					iter = examples.begin();
					std::random_shuffle(iter, iter_end);
					algo.OpenEpoch(false);
					for (; iter != iter_end; ++iter)
						algo.Submit(*iter);
					algo.CloseEpoch();
				}

				void labelize() {
					typename GNG_T::Edges::iterator edge_iter, edge_end;
					typename GNG_T::Nodes::iterator node_iter, node_end;
					typename GNG_T::Node *n;
					typename GNG_T::Edge *e;
					double length_max;

					length_max = param_gngt_length_max()
							* param_gngt_length_max();

					for (edge_iter = algo.edges.begin(), edge_end
							= algo.edges.end(); edge_iter != edge_end; ++edge_iter) {
						e = *edge_iter;
						labelizer.SetValidity(e, Input::d2(e->n1->value.w,
								e->n2->value.w) < length_max);
					}

					for (node_iter = algo.nodes.begin(), node_end
							= algo.nodes.end(); node_iter != node_end; ++node_iter) {
						n = *node_iter;
						labelizer.SetValidity(n, n->value.e / n->value.n
								< param_gngt_variance_max());
					}
					labelizer.Process(algo);
				}

				void draw(ImageRGB24*& result, ImageRGB24& original) {
					typename vq::Labelizer<GNG_T>::Labeling::iterator
							label_iter, label_end;
					typename vq::Labelizer<GNG_T>::ConnectedComponent
							* component;
					typename vq::Labelizer<GNG_T>::ConnectedComponent::Edges::iterator
							edge_iter, edge_end;
					int label;

					mirage::img::Line<ImageRGB24> line;
					mirage::img::Line<ImageRGB24>::pixel_type lpix, lpix_end;
					mirage::img::Coordinate A, B;
					typename GNG_T::Node *n1, *n2;
					mirage::colorspace::RGB_24 paint;
					mirage::img::Coordinate img_origin, img_size;

					result->resize(original._dimension);
					*result = original;
					mirage::SubFrame<ImageRGB24> dot(*result, pen, pen); // silly pen args

					line << (*result);
					for (label_iter = labelizer.labeling.begin(), label_end
							= labelizer.labeling.end(); label_iter != label_end; ++label_iter) {

						label = label_iter->first;
						component = label_iter->second;

						if (colors.count(label) == 0)
							colors[label] = getRandomColor();
						paint = colors[label];

						for (edge_iter = component->edges.begin(), edge_end
								= component->edges.end(); edge_iter != edge_end; ++edge_iter) {
							n1 = (*edge_iter)->n1;
							n2 = (*edge_iter)->n2;

							A((int) (n1->value.w[0] + .5),
									(int) (n1->value.w[1] + .5));
							B((int) (n2->value.w[0] + .5),
									(int) (n2->value.w[1] + .5));
							if ((A[0] != B[0]) || (A[1] != B[1])) {
								line(A, B, false, false);
								for (lpix = line.begin(), lpix_end = line.end(); lpix
										!= lpix_end; ++lpix) {
									dot.resize(!lpix - half_pen, pen);
									dot = paint;
								}
							}
						}
					}
				}

				void draw(ImageRGB24*& result, ImageBool& contour) {
					typename vq::Labelizer<GNG_T>::Labeling::iterator
							label_iter, label_end;
					typename vq::Labelizer<GNG_T>::ConnectedComponent
							* component;
					typename vq::Labelizer<GNG_T>::ConnectedComponent::Edges::iterator
							edge_iter, edge_end;
					int label;

					mirage::img::Line<ImageRGB24> line;
					mirage::img::Line<ImageRGB24>::pixel_type lpix, lpix_end;
					mirage::img::Coordinate A, B;
					typename GNG_T::Node *n1, *n2;
					mirage::colorspace::RGB_24 paint;
					mirage::img::Coordinate img_origin, img_size;

					result->resize(contour._dimension);
					*result = mirage::colorspace::RGB_24(0, 0, 0);

					mirage::SubFrame<ImageRGB24> dot(*result, pen, pen); // silly pen args

					line << (*result);
					for (label_iter = labelizer.labeling.begin(), label_end
							= labelizer.labeling.end(); label_iter != label_end; ++label_iter) {

						label = label_iter->first;
						component = label_iter->second;

						if (colors.count(label) == 0)
							colors[label] = getRandomColor();
						paint = colors[label];

						for (edge_iter = component->edges.begin(), edge_end
								= component->edges.end(); edge_iter != edge_end; ++edge_iter) {
							n1 = (*edge_iter)->n1;
							n2 = (*edge_iter)->n2;

							A((int) (n1->value.w[0] + .5),
									(int) (n1->value.w[1] + .5));
							B((int) (n2->value.w[0] + .5),
									(int) (n2->value.w[1] + .5));
							if ((A[0] != B[0]) || (A[1] != B[1])) {
								line(A, B, false, false);
								for (lpix = line.begin(), lpix_end = line.end(); lpix
										!= lpix_end; ++lpix) {
									dot.resize(!lpix - half_pen, pen);
									dot = paint;
								}
							}
						}
					}
				}

			public:
				/** initialize drawing colors */
				GNGT_Draw() {
					pen[0] = pen[1] = param_pen_thickness();

					half_pen = pen / 2;

					colors[0]._red = 255;
					colors[0]._green = 95;
					colors[0]._blue = 95;

					colors[1]._red = 0;
					colors[1]._green = 0;
					colors[1]._blue = 255;

					colors[2]._red = 0;
					colors[2]._green = 255;
					colors[2]._blue = 0;

					colors[3]._red = 0;
					colors[3]._green = 255;
					colors[3]._blue = 255;

					colors[4]._red = 190;
					colors[4]._green = 190;
					colors[4]._blue = 255;

					colors[5]._red = 255;
					colors[5]._green = 0;
					colors[5]._blue = 255;

					colors[6]._red = 255;
					colors[6]._green = 190;
					colors[6]._blue = 0;

					colors[7]._red = 190;
					colors[7]._green = 255;
					colors[7]._blue = 0;

					colors[8]._red = 255;
					colors[8]._green = 190;
					colors[8]._blue = 190;
				}

				/**
				 * @param contour the image to post
				 * @param outFrame the original image with the gngt graph drawn over it
				 * @param fromQueue image received from ImageFeeder via a bypass queue - it contains the original image
				 */
				void operator()(ImageBool* contour, ImageRGB24*& outFrame,
						ImageRGB24& fromQueue) {
					feedGNGT(*contour);
					labelize();
					draw(outFrame, fromQueue);
				}

				/**
				 * @param contour the image to post
				 * @param image containing the gngt graph
				 */
				void operator()(ImageBool* contour, ImageRGB24*& outFrame) {
					feedGNGT(*contour);
					labelize();
					draw(outFrame, *contour);
				}

		};
};

#endif /* TRACK_H_ */
