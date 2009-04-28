#include "track.h"

using namespace track;

#define KWD_BACKGROUND_SAMPLE "BACKGROUND_SAMPLES"
#define KWD_BACKGROUND_COEF   "BACKGROUND_COEF"
#define KWD_INTERFRAME_DELAY  "FRAME_DELAY"
#define KWD_BUFFER_IN         "BUFFER_IN"
#define KWD_BUFFER_OUT        "BUFFER_OUT"
#define KWD_QUALITY_OUT       "JPEG_QUALITY_OUT"
#define KWD_BUFFER_DETECTION  "BUFFER_DETECTION"
#define KWD_QUALITY_DETECTION "JPEG_QUALITY_DETECTION"
#define KWD_DISPLAY_DETECTION "DISPLAY_DETECTION"
#define KWD_COS               "COLOR_MIN_COSINE"
#define KWD_NORM              "COLOR_MIN_SQUARE_NORM"
#define KWD_NORM_MARGIN       "COLOR_MIN_SQUARE_NORM_MARGIN"
#define KWD_TARGET            "GNGT_TARGET"
#define KWD_AGE               "GNGT_EDGE_AGE"
#define KWD_LEARNING_RATES    "GNGT_LEARNING_RATES"
#define KWD_VARIANCE          "GNGT_MAX_VARIANCE"
#define KWD_LENGTH            "GNGT_MAX_LENGTH"
#define KWD_EPOCHS            "GNGT_EPOCHS_PER_FRAME"
#define KWD_PEN_THICKNESS     "PEN_THICKNESS"
#define KWD_MARGINS           "MARGINS"
#define KWD_MORPH             "MORPHOMATH"

/*-----------------------------------------------------------------------------------------*/
/* define the global vars in track::                                                       */
/*-----------------------------------------------------------------------------------------*/
int track::param_inter_frame_delay;
std::string track::param_buffer_in_hostname;
std::string track::param_buffer_in_resource;
std::string track::param_buffer_in_entity;
int track::param_buffer_in_port;
bool track::param_buffer_in_jpg;
std::string track::param_buffer_out_hostname;
std::string track::param_buffer_out_resource;
std::string track::param_buffer_out_entity;
int track::param_buffer_out_port;
bool track::param_buffer_out_jpg;
int track::param_out_jpeg_quality;
std::string track::param_buffer_detection_hostname;
std::string track::param_buffer_detection_resource;
std::string track::param_buffer_detection_entity;
int track::param_buffer_detection_port;
bool track::param_buffer_detection_jpg;
int track::param_detection_jpeg_quality;
bool track::param_display_detection;
double track::param_color_min_norm2;
double track::param_color_min_norm2_margin;
double track::param_color_min_cosine;
int track::param_nb_background_samples;
double track::param_background_coef;
double track::param_gngt_target;
double track::param_gngt_first_learning_rate;
double track::param_gngt_second_learning_rate;
int track::param_gngt_edge_age_max;
double track::param_gngt_variance_max;
double track::param_gngt_length_max;
int track::param_nb_epochs_per_frame;
int track::param_pen_thickness;
int track::param_left_margin;
int track::param_right_margin;
int track::param_top_margin;
int track::param_bottom_margin;
bool track::param_morph;
int track::param_morph_radius;

ServerConnection* track::serverConn;

/*-----------------------------------------------------------------------------------------*/
/* set default parameters                                                                  */
/*-----------------------------------------------------------------------------------------*/
void track::ParameterParser::loadDefaultParameters() {
	param_nb_background_samples = 10;
	param_inter_frame_delay = 10;
	param_buffer_in_hostname = "localhost";
	param_buffer_in_resource = "JPEG-RESOURCE";
	param_buffer_in_entity = "video-buf";
	param_buffer_in_port = 20000;
	param_buffer_in_jpg = true;
	param_buffer_out_hostname = "localhost";
	param_buffer_out_resource = "JPEG-RESOURCE";
	param_buffer_out_entity = "display-buf1";
	param_buffer_out_port = 20000;
	param_buffer_out_jpg = true;
	param_out_jpeg_quality = 70;
	param_buffer_detection_hostname = "localhost";
	param_buffer_detection_resource = "JPEG-RESOURCE";
	param_buffer_detection_entity = "display-buf2";
	param_buffer_detection_port = 20000;
	param_buffer_detection_jpg = true;
	param_detection_jpeg_quality = 70;
	param_display_detection = false;
	param_color_min_norm2 = 7000;
	param_color_min_norm2_margin = 100;
	param_color_min_cosine = .997;
	param_background_coef = .0005;
	param_gngt_target = 30;
	param_gngt_first_learning_rate = .05;
	param_gngt_second_learning_rate = .005;
	param_gngt_edge_age_max = 20;
	param_gngt_variance_max = 50;
	param_gngt_length_max = 20;
	param_nb_epochs_per_frame = 10;
	param_pen_thickness = 3;
	param_top_margin = 5;
	param_bottom_margin = 5;
	param_left_margin = 5;
	param_right_margin = 5;
	param_morph = false;
	param_morph_radius = 3;
}
track::ParameterParser::ParameterParser(void) {
	loadDefaultParameters();
	serverConn = new ServerConnection();
}

/*-----------------------------------------------------------------------------------------*/
/* get the parameters from the input file                                                  */
/*-----------------------------------------------------------------------------------------*/
track::ParameterParser::ParameterParser(const std::string& filename) {
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
			if (kwd == KWD_BACKGROUND_SAMPLE)
				file >> param_nb_background_samples;
			else if (kwd == KWD_BACKGROUND_COEF)
				file >> param_background_coef;
			else if (kwd == KWD_INTERFRAME_DELAY)
				file >> param_inter_frame_delay;
			else if (kwd == KWD_BUFFER_IN) {
				file >> param_buffer_in_hostname >> param_buffer_in_port
						>> param_buffer_in_resource >> param_buffer_in_entity
						>> mode;
				param_buffer_in_jpg = true || mode == "jpg" || mode == "jpeg"
						|| mode == "JPG" || mode == "JPEG";

			} else if (kwd == KWD_BUFFER_OUT) {
				file >> param_buffer_out_hostname >> param_buffer_out_port
						>> param_buffer_out_resource >> param_buffer_out_entity
						>> mode;
				param_buffer_out_jpg = true || mode == "jpg" || mode == "jpeg"
						|| mode == "JPG" || mode == "JPEG";

			} else if (kwd == KWD_QUALITY_OUT)
				file >> param_out_jpeg_quality;
			else if (kwd == KWD_BUFFER_DETECTION) {
				file >> param_buffer_detection_hostname
						>> param_buffer_detection_port
						>> param_buffer_detection_resource
						>> param_buffer_detection_entity >> mode;
				param_buffer_detection_jpg = true || mode == "jpg" || mode
						== "jpeg" || mode == "JPG" || mode == "JPEG";
			} else if (kwd == KWD_QUALITY_DETECTION)
				file >> param_detection_jpeg_quality;
			else if (kwd == KWD_DISPLAY_DETECTION)
				param_display_detection = true;
			else if (kwd == KWD_COS)
				file >> param_color_min_cosine;
			else if (kwd == KWD_NORM_MARGIN)
				file >> param_color_min_norm2_margin;
			else if (kwd == KWD_NORM)
				file >> param_color_min_norm2;
			else if (kwd == KWD_TARGET)
				file >> param_gngt_target;
			else if (kwd == KWD_AGE)
				file >> param_gngt_edge_age_max;
			else if (kwd == KWD_LEARNING_RATES)
				file >> param_gngt_first_learning_rate
						>> param_gngt_second_learning_rate;
			else if (kwd == KWD_VARIANCE)
				file >> param_gngt_variance_max;
			else if (kwd == KWD_LENGTH)
				file >> param_gngt_length_max;
			else if (kwd == KWD_EPOCHS)
				file >> param_nb_epochs_per_frame;
			else if (kwd == KWD_PEN_THICKNESS)
				file >> param_pen_thickness;
			else if (kwd == KWD_MARGINS)
				file >> param_left_margin >> param_right_margin
						>> param_top_margin >> param_bottom_margin;
			else if (kwd == KWD_MORPH) {
				param_morph = true;
				file >> param_morph_radius;
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

	serverConn = new ServerConnection();
}

/*-----------------------------------------------------------------------------------------*/
/* save current parameters to file                                                         */
/*-----------------------------------------------------------------------------------------*/
void track::ParameterParser::saveParameters(const std::string& filename) {
	std::ofstream file;

	file.open(filename.c_str());
	if (!file) {
		std::cerr << "Cannot open file \"" << filename
				<< "\" for saving configuration." << std::endl;
		::exit(1);
	}
	file << "######################################" << std::endl
			<< "#                                    #" << std::endl
			<< "# Parameters for gngt-video tracking #" << std::endl
			<< "#                                    #" << std::endl
			<< "######################################" << std::endl
			<< std::endl
			<< "# Margins of the computation area from sides : left, right, top and bottom."
			<< std::endl << KWD_MARGINS << ' ' << param_left_margin << ' '
			<< param_right_margin << ' ' << param_top_margin << ' '
			<< param_bottom_margin << std::endl << std::endl
			<< "# Number of sample frames for background initialization."
			<< std::endl << KWD_BACKGROUND_SAMPLE << ' '
			<< param_nb_background_samples << std::endl
			<< "# Background updating low-pass filter coefficient."
			<< std::endl << KWD_BACKGROUND_COEF << ' ' << param_background_coef
			<< std::endl << std::endl << "# Interframe grabbing delay (ms)."
			<< std::endl << KWD_INTERFRAME_DELAY << ' '
			<< param_inter_frame_delay << std::endl << std::endl
			<< "# Input and output. Order is host, port, resource, entity, mode (JPEG | IMG)."
			<< std::endl << KWD_BUFFER_IN << ' ' << param_buffer_in_hostname
			<< ' ' << param_buffer_in_port << ' ' << param_buffer_in_resource
			<< ' ' << param_buffer_in_entity << ' ' << Mode(param_buffer_in_jpg) << std::endl << KWD_BUFFER_OUT << ' '
			<< param_buffer_out_hostname << ' ' << param_buffer_out_port << ' '
			<< param_buffer_out_resource << ' ' << param_buffer_out_entity
			<< ' ' << Mode(param_buffer_out_jpg) << std::endl << std::endl
			<< "# Idf output is compressed, this is its Jpeg quality (1--100)."
			<< std::endl << KWD_QUALITY_OUT << ' ' << param_out_jpeg_quality
			<< std::endl << std::endl
			<< "# Comment out the following to disable the sending"
			<< std::endl << "# of object detection result to server."
			<< std::endl;
	if (!param_display_detection)
		file << "# ";
	file << KWD_DISPLAY_DETECTION << std::endl
			<< "# If enabled, this describe the sending of detection to the server."
			<< std::endl << KWD_BUFFER_DETECTION << ' '
			<< param_buffer_detection_hostname << ' '
			<< param_buffer_detection_port << ' '
			<< param_buffer_detection_resource << ' '
			<< param_buffer_detection_entity << ' ' << Mode(
			param_buffer_detection_jpg) << std::endl << KWD_QUALITY_DETECTION
			<< ' ' << param_detection_jpeg_quality << std::endl << std::endl
			<< "# Color comparison : minimal square norm, its tolerance,"
			<< std::endl << "# and the minimal cosine." << std::endl
			<< KWD_NORM << ' ' << param_color_min_norm2 << std::endl
			<< KWD_NORM_MARGIN << ' ' << param_color_min_norm2_margin
			<< std::endl << KWD_COS << ' ' << param_color_min_cosine
			<< std::endl << std::endl
			<< "# GNG-T Parameters : target, maximum edge age," << std::endl
			<< "# winner and second learning rates." << std::endl << KWD_TARGET
			<< ' ' << param_gngt_target << std::endl << KWD_AGE << ' '
			<< param_gngt_edge_age_max << std::endl << KWD_LEARNING_RATES
			<< ' ' << param_gngt_first_learning_rate << ' '
			<< param_gngt_second_learning_rate << std::endl
			<< "# Number of GNG-T epochs per frame (a supplementary one will be actually used)."
			<< std::endl << KWD_EPOCHS << ' ' << param_nb_epochs_per_frame
			<< std::endl << "# Maximal variance for a non noisy node."
			<< std::endl << KWD_VARIANCE << ' ' << param_gngt_variance_max
			<< std::endl << "# Maximal edge length inside a polygon."
			<< std::endl << KWD_LENGTH << ' ' << param_gngt_length_max
			<< std::endl << std::endl << "# Drawing paremeters" << std::endl
			<< KWD_PEN_THICKNESS << ' ' << param_pen_thickness << std::endl
			<< std::endl
			<< "# Uncomment the following line to enable morpho-matematical"
			<< std::endl << "# detection cleaning. Number is the mask radius."
			<< std::endl;
	if (!param_morph)
		file << "# ";
	file << KWD_MORPH << ' ' << param_morph_radius << std::endl << std::endl
			<< std::endl;

	file.close();
}

/*-----------------------------------------------------------------------------------------*/
/* connect to server using the parameters in the track namespace                           */
/*-----------------------------------------------------------------------------------------*/
track::ServerConnection::ServerConnection() {
	this->jpgClientDetection = NULL;
	this->jpgClientIn = NULL;
	this->jpgClientOut = NULL;
	this->rawClientDetection = NULL;
	this->rawClientIn = NULL;
	this->rawClientOut = NULL;

	std::cout << "  --> " << param_buffer_in_entity << '('
			<< param_buffer_in_resource << ")@" << param_buffer_in_hostname
			<< ':' << param_buffer_in_port << std::endl;
	if (param_buffer_in_jpg) {
		jpgClientIn = new bkbd::Repository<bkbd::JPEG>::Client(
				param_buffer_in_resource);
		if (!jpgClientIn->Connect(param_buffer_in_hostname,
				param_buffer_in_port)) {
			std::cerr << "Cannot connect to \"" << param_buffer_in_resource
					<< "\"@" << param_buffer_in_hostname << ':'
					<< param_buffer_in_port << std::endl;
			::exit(1);
		}
	} else {
		rawClientIn = new bkbd::Repository<bkbd::Image>::Client(
				param_buffer_in_resource);
		if (!rawClientIn->Connect(param_buffer_in_hostname,
				param_buffer_in_port)) {
			std::cerr << "Cannot connect to \"" << param_buffer_in_resource
					<< "\"@" << param_buffer_in_hostname << ':'
					<< param_buffer_in_port << std::endl;
			::exit(1);
		}
	}

	std::cout << "  <-- " << param_buffer_out_entity << '('
			<< param_buffer_out_resource << ")@" << param_buffer_out_hostname
			<< ':' << param_buffer_out_port << std::endl;
	if (param_buffer_out_jpg) {
		jpgClientOut = new bkbd::Repository<bkbd::JPEG>::Client(
				param_buffer_out_resource);
		if (!jpgClientOut->Connect(param_buffer_out_hostname,
				param_buffer_out_port)) {
			std::cerr << "Cannot connect to \"" << param_buffer_out_resource
					<< "\"@" << param_buffer_out_hostname << ':'
					<< param_buffer_out_port << std::endl;
			::exit(1);
		}
	} else {
		rawClientOut = new bkbd::Repository<bkbd::Image>::Client(
				param_buffer_out_resource);
		if (!rawClientOut->Connect(param_buffer_out_hostname,
				param_buffer_out_port)) {
			std::cerr << "Cannot connect to \"" << param_buffer_out_resource
					<< "\"@" << param_buffer_out_hostname << ':'
					<< param_buffer_out_port << std::endl;
			::exit(1);
		}
	}

	if (param_display_detection) {
		std::cout << "  <-- " << param_buffer_detection_entity << '('
				<< param_buffer_detection_resource << ")@"
				<< param_buffer_detection_hostname << ':'
				<< param_buffer_detection_port << std::endl;
		if (param_buffer_detection_jpg) {
			jpgClientDetection = new bkbd::Repository<bkbd::JPEG>::Client(
					param_buffer_detection_resource);
			if (!jpgClientDetection->Connect(param_buffer_detection_hostname,
					param_buffer_detection_port)) {
				std::cerr << "Cannot connect to \""
						<< param_buffer_detection_resource << "\"@"
						<< param_buffer_detection_hostname << ':'
						<< param_buffer_detection_port << std::endl;
				::exit(1);
			}
		} else {
			rawClientDetection = new bkbd::Repository<bkbd::Image>::Client(
					param_buffer_detection_resource);
			if (!rawClientDetection->Connect(param_buffer_detection_hostname,
					param_buffer_detection_port)) {
				std::cerr << "Cannot connect to \""
						<< param_buffer_detection_resource << "\"@"
						<< param_buffer_detection_hostname << ':'
						<< param_buffer_detection_port << std::endl;
				::exit(1);
			}
		}
	} else
		std::cout << "  <-- detection display disabled." << std::endl;

}

track::ServerConnection::~ServerConnection() {
	//TODO close connection
}

/*-----------------------------------------------------------------------------------------*/
/* retrieve an image from the server                                                       */
/*-----------------------------------------------------------------------------------------*/
void track::ServerConnection::getImage(bkbd::Image* frame) {
	static bkbd::JPEG in_jpeg_buf;
	static struct timeval timestamp;
	static jpeg::Decompress decompressor;
	int w, h, d;

	frame->setAllocationPolicy(bkbd::Image::AllocateAutomatic);
	if (param_buffer_in_jpg) {
		jpgClientIn->get(param_buffer_in_entity, in_jpeg_buf, timestamp);
		decompressor.setInputStream(in_jpeg_buf.ijpeg());
		decompressor.Flush();
		decompressor.readHeader(w, h, d);
		frame->resize(w, h, (bkbd::Image::Depth) d);
		decompressor.readImage(frame->data);
	} else
		rawClientIn->get(param_buffer_in_entity, *frame, timestamp);

}

/*-----------------------------------------------------------------------------------------*/
/* post image to the server                                                                */
/*-----------------------------------------------------------------------------------------*/
void track::ServerConnection::postImg(ImageRGB24& outFrame) {
	static bkbd::Image dest;
	static jpeg::Compress compressor;
	static bkbd::JPEG out_jpeg_buf;

	dest.setAllocationPolicy(bkbd::Image::AllocateNone,
			(unsigned char*) (&(*outFrame.begin())));
	dest.resize(outFrame._dimension[0], outFrame._dimension[1],
			bkbd::Image::RGB24);

	if (param_buffer_out_jpg) {
		compressor.setOutputStream(out_jpeg_buf.ojpeg());
		compressor.writeImage(dest.width, dest.height, (int) (dest.depth),
				dest.data, param_out_jpeg_quality);
		jpgClientOut->post(param_buffer_out_entity, out_jpeg_buf);
	} else
		rawClientOut->post(param_buffer_out_entity, dest);
}

/*-----------------------------------------------------------------------------------------*/
/* retrieve images from the server                                                         */
/*-----------------------------------------------------------------------------------------*/
void track::ImageFeeder::operator ()(void* dummy, bkbd::Image*& outFrame) {
	track::serverConn->getImage(outFrame);
}

/*-----------------------------------------------------------------------------------------*/
/* post images to the server                                                               */
/*-----------------------------------------------------------------------------------------*/
void track::ImagePoster::operator ()(ImageRGB24* image, void* dummy) {
	track::serverConn->postImg(*image);
}

/*-----------------------------------------------------------------------------------------*/
/* obtain grayscale image from bkbd::Image                                                 */
/*-----------------------------------------------------------------------------------------*/
void track::Convert2RGB::operator ()(bkbd::Image* inFrame, ImageRGB24*& outFrame) {
	outFrame->resize(mirage::img::Coordinate(inFrame->width, inFrame->height),
			(mirage::colorspace::RGB_24*) (inFrame->data));
}

void track::Convert2Gray::operator ()(ImageRGB24* inFrame, ImageGRAY8*& outFrame) {
	outFrame->resize(inFrame->_dimension);
	mirage::algo::UnaryOp<ImageRGB24, ImageGRAY8,
			mirage::colorspace::RGBToGray<ImageRGB24::value_type,
					ImageGray8::value_type> >(*inFrame, *outFrame);
}

/*-----------------------------------------------------------------------------------------*/
/* foreground detection                                                                    */
/*-----------------------------------------------------------------------------------------*/
void track::ForeGround::operator()(ImageGRAY8* inFrame, ImageBool*& outFrame) {
	if (firstImg) {
		algo.initForegroundDetection(*inFrame);
		firstImg = false;
	}
	algo.SigmaDeltaModified(*inFrame);
	algo.getMask(*outFrame);
}

/*-----------------------------------------------------------------------------------------*/
/* eliminate noise from the image                                                          */
/*-----------------------------------------------------------------------------------------*/
void track::MorphoMath::operator ()(ImageBool* inFrame, ImageBool*& outFrame) {
	ImageBool::pixel_type pix1, pix2, pix_end;
	ImageBool::point_type pos, offset;
	bool n, s, e, w;

	outFrame->resize(inFrame->_dimension);
	if (!param_morph) {
		for (pix1 = inFrame->begin(), pix_end = inFrame->end(), pix2
				= outFrame->begin(); pix1 != pix_end; ++pix1, ++pix2) {
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
		mirage::morph::Format<ImageBool, ImageBool, 0>::Opening(*inFrame,
				element, *outFrame);
	}
}

/*-----------------------------------------------------------------------------------------*/
/* get contour from the foreground image                                                   */
/*-----------------------------------------------------------------------------------------*/
void track::Contour::operator ()(ImageBool* inFrame, ImageBool*& outFrame) {
	ImageBool::pixel_type pix1, pix2, pix_end;
	ImageBool::point_type pos, dimension, offset;

	outFrame->resize(inFrame->_dimension);
	for (pix1 = inFrame->begin(), pix_end = inFrame->end(), pix2
			= outFrame->begin(); pix1 != pix_end; ++pix1, ++pix2) {
		if (*pix1) {
			pos = !pix2 + offset(1, 0);
			*pix2 = !((*inFrame)(pos));

			if (!(*pix2)) {
				pos = !pix2 + offset(0, 1);
				*pix2 = !((*inFrame)(pos));

				if (!(*pix2)) {

					pos = !pix2 + offset(0, -1);
					*pix2 = !((*inFrame)(pos));

					if (!(*pix2)) {
						pos = !pix2 + offset(-1, 0);
						*pix2 = !((*outFrame)(pos));
					}
				}
			}
		} else {
			*pix2 = false;
		}
	}
}

/*-----------------------------------------------------------------------------------------*/
/* apply gngt on the contour                                                               */
/*-----------------------------------------------------------------------------------------*/
void track::GNGT::feedGNGT(ImageBool& contour) {
	ImageBool::pixel_type pix, pix_end;
	std::vector<Input>::iterator iter, iter_end;
	Input example;
	int epoch;

	// Getting examples
	examples.clear();
	for (pix = contour.begin(), pix_end = contour.end(); pix != pix_end; ++pix)
		if (*pix) {
			example[0] = (!pix)[0];
			example[1] = (!pix)[1];
			examples.push_back(example);
		}

	for (epoch = 0; epoch < param_nb_epochs_per_frame; ++epoch) {
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

void track::GNGT::labelize(LABELIZER& labelizer) {
	GNG_T::Edges::iterator edge_iter, edge_end;
	GNG_T::Nodes::iterator node_iter, node_end;
	GNG_T::Node *n;
	GNG_T::Edge *e;
	double length_max;

	length_max = param_gngt_length_max * param_gngt_length_max;

	for (edge_iter = algo.edges.begin(), edge_end = algo.edges.end(); edge_iter
			!= edge_end; ++edge_iter) {
		e = *edge_iter;
		labelizer.SetValidity(e, Input::d2(e->n1->value.w, e->n2->value.w)
				< length_max);
	}

	for (node_iter = algo.nodes.begin(), node_end = algo.nodes.end(); node_iter
			!= node_end; ++node_iter) {
		n = *node_iter;
		labelizer.SetValidity(n, n->value.e / n->value.n
				< param_gngt_variance_max);
	}
	labelizer.Process(algo);
}

/*-----------------------------------------------------------------------------------------*/
/* draw the graph over the original image                                                  */
/*-----------------------------------------------------------------------------------------*/
mirage::colorspace::RGB_24 track::Draw::getRandomColor() {
	mirage::colorspace::RGB_24 res;

	res._red = (mirage::colorspace::RGB_24::value_type) (256.0 * (rand()
			/ (RAND_MAX + 1.0)));
	res._green = (mirage::colorspace::RGB_24::value_type) (256.0 * (rand()
			/ (RAND_MAX + 1.0)));
	res._blue = (mirage::colorspace::RGB_24::value_type) (256.0 * (rand()
			/ (RAND_MAX + 1.0)));

	return res;
}

/* initialize colors */
track::Draw::Draw() {
	pen[0] = pen[1] = param_pen_thickness;
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

void track::Draw::operator ()(LABELIZER* labelizer, ImageRGB24*& result, ImageRGB24& original) {
	vq::Labelizer<GNG_T>::Labeling::iterator label_iter, label_end;
	vq::Labelizer<GNG_T>::ConnectedComponent* component;
	vq::Labelizer<GNG_T>::ConnectedComponent::Edges::iterator edge_iter,
			edge_end;
	int label;
	mirage::img::Line<ImageRGB24> line;
	mirage::img::Line<ImageRGB24>::pixel_type lpix, lpix_end;
	mirage::SubFrame<ImageRGB24> dot(*result, pen, pen);
	mirage::img::Coordinate A, B;
	GNG_T::Node *n1, *n2;
	mirage::colorspace::RGB_24 paint;

	mirage::SubFrame<ImageRGB24> source(original,
			mirage::img::Coordinate(0, 0), mirage::img::Coordinate(0, 0));

	result->resize(original._dimension);

	std::cout << source.size() << " " << result->size() << " "
			<< original.size() << "\n";

	mirage::algo::UnaryOp<mirage::SubFrame<ImageRGB24>, ImageRGB24,
			mirage::algo::Affectation<mirage::SubFrame<ImageRGB24>::value_type,
					ImageRGB24::value_type> >(source, *result);

	line << (*result);
	for (label_iter = labelizer->labeling.begin(), label_end
			= labelizer->labeling.end(); label_iter != label_end; ++label_iter) {

		label = label_iter->first;
		component = label_iter->second;

		if (colors.count(label) == 0)
			colors[label] = getRandomColor();
		paint = colors[label];

		for (edge_iter = component->edges.begin(), edge_end
				= component->edges.end(); edge_iter != edge_end; ++edge_iter) {
			n1 = (*edge_iter)->n1;
			n2 = (*edge_iter)->n2;
			A((int) (n1->value.w[0] + .5), (int) (n1->value.w[1] + .5));
			B((int) (n2->value.w[0] + .5), (int) (n2->value.w[1] + .5));
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
