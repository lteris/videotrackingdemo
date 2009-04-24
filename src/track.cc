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
/* set default parameters                                                                  */
/*-----------------------------------------------------------------------------------------*/
track::ParameterParser::ParameterParser(void) {
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
	SetDefaultParameters();

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
			<< ' ' << param_buffer_in_entity << ' '
			<< Mode(param_buffer_in_jpg) << std::endl << KWD_BUFFER_OUT << ' '
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

	src.setAllocationPolicy(bkbd::Image::AllocateAutomatic);
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
void track::ServerConnection::postImg(const ImageRGB24& outFrame) {
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
		raw_client_out->post(param_buffer_out_entity, dest);
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
	} else {
		algo.SigmaDeltaModified(*inFrame);
		algo.getMask(*outFrame);
	}
}

/*-----------------------------------------------------------------------------------------*/
/* eliminate noise from the image                                                          */
/*-----------------------------------------------------------------------------------------*/
void track::MorphoMath::operator ()(ImageBool* inFrame, ImageBool*& outFrame) {
	ImageBool::pixel_type pix1, pix2, pix_end;
	ImageBool::point_type pos, offset;
	bool n, s, e, w;

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
			} else
				*pix2 = false;
		}
	} else
		mirage::morph::Format<ImageBool, ImageBool, 0>::Opening(src, element,
				res);
}
