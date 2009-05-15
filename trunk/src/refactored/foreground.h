/*   This file is part of gngtlib
 *
 *   Copyright (C) 2009,  Supelec
 *
 *   Author : Herve Frezza-Buet
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
#ifndef mirageFOREGROUND_H
#define mirageFOREGROUND_H

#include <mirage.h>
#include <math.h>

//Mirage typedefs
typedef mirage::img::Coding<mirage::colorspace::RGB_24>::Frame ImageRGB24;
typedef mirage::img::Coding<mirage::colorspace::GRAY_8>::Frame ImageGray8;
typedef mirage::img::Coding<bool>::Frame                       ImageBool;
typedef mirage::img::Coding<mirage::morph::MaskValue>::Frame   Element;

namespace mirage {

  namespace img {

    /**
     * @short Foreground estimation class.
     * This class computes an estimation of the foreground objects in a series
     * of images
     * 
     */
    template <typename PARAM, int MORPH=0>
    class ForegroundDetection {
    public:
    
    /**
     * Constructor for the foreground detection object.
     * @param imgExample The first image of the image series to intialize the estimation.
     */
    void initForegroundDetection(ImageGray8& imgExample) {
      imgSD_M.resize(imgExample._dimension);
      imgSD_D.resize(imgExample._dimension);
      imgSD_V.resize(imgExample._dimension);
      imgSD_U_M.resize(imgExample._dimension); 
      imgSD_U_STD.resize(imgExample._dimension);
      CopyImage(imgExample, imgSD_M);
      CopyImage(imgExample, imgSD_D);
      CopyImage(imgExample, imgSD_V);
      MakeNull(imgSD_U_M, imgSD_U_M);
      MakeNull(imgSD_U_STD, imgSD_U_STD);
      InitImage(imgSD_V, imgSD_V);
      
      imgBoolOut.resize(imgExample._dimension);
      imgBoolClean.resize(imgExample._dimension);
      imgBoolMorph.resize(imgExample._dimension);
      //Initialise the morpho disks if we want morpho      
      if (MORPH){
	mirage::morph::Mask::Disk(elementDil,PARAM::morphDil());
	mirage::morph::Mask::Disk(elementEro,PARAM::morphEro());
      }
    }
    
    
    /**
     * Calculates the image mask representing the foreground estimation using a modified
     * Sigma-Delta algorithm.
     * @param imgIn The next image in the series.
     */
    void SigmaDeltaModified(ImageGray8& imgIn) {
      ImageGray8::pixel_type        pixM, pixD, pixV, pixU_M, pixU_STD;
      ImageGray8::const_pixel_type  cpix, cpix_end;
      ImageBool::pixel_type         pix, pix_end;
      int                           temp, vMax, vMin;
      double                        dTemp, k;
      
      vMax = PARAM::stdMax();
      vMin = PARAM::stdMin();
      k = PARAM::k();
      
      for(pix = imgBoolOut.begin(), pix_end = imgBoolOut.end(),
	    cpix = imgIn.const_begin(), cpix_end = imgIn.const_end(),
	    pixM = imgSD_M.begin(), pixV = imgSD_V.begin(),
	    pixD = imgSD_D.begin(), pixU_M = imgSD_U_M.begin(),pixU_STD = imgSD_U_STD.begin();
	  pix != pix_end && cpix != cpix_end;
	  ++pix,++cpix,++pixM,++pixV,++pixD,++pixU_M, ++pixU_STD) {
	
	*pixD = abs(*pixM - *cpix);
	dTemp = k * (double) (*pixV);
	
	if ((double) *pixD > dTemp) {
	  
	  if (*pixU_M == 1) {
	    temp = *cpix - *pixM;
	    if (temp > 0) {
	      *pixM = *pixM + 1;
	    } else if (temp > 0) {
	      *pixM = *pixM - 1;
	    }
	    *pixU_M = 0;
	  } else if (*pixU_M == 0) {
	    temp = 4 * (*pixV);
	    if (*pixD > temp) {
	      *pixU_M = 2;
	    } else {
	      *pixU_M = 2;
	    }
	  } else if (*pixU_M > 1) {
	    *pixU_M = *pixU_M - 1;
	  }
	  *pix = 1;
	  
	} else {
	  
	  temp = *cpix - *pixM;
	  if (temp > 0) {
	    *pixM = *pixM + 1;
	  } else if (temp > 0) {
	    *pixM = *pixM - 1;
	  }
	  *pixU_M = 0;
	  
	  if (*pixD > *pixV) {

	    *pixV = *pixV + 1;

	  } else if(*pixD < *pixV) {

	    *pixU_STD = *pixU_STD - 1;
	    if ((*pixU_STD == 0) && (*pixV > vMin)) {
	      *pixV = *pixV - 1;
	      if (*pixV > vMax) {
		*pixU_STD = 4;
	      } else {
		*pixU_STD = 6;
	      }
	    }
 
	  }
	  
	  *pix = 0;
	}
      }
      
      ShadowRemove(imgIn, imgBoolOut);
      CleanImage(imgBoolOut, imgBoolClean);

      if (MORPH) {
	//Do morphological operations
	mirage::morph::Format<ImageBool,ImageBool,1>::Erosion(imgBoolOut,elementEro,imgBoolMorph);
	mirage::morph::Format<ImageBool,ImageBool,1>::Dilatation(imgBoolOut,elementDil,imgBoolMorph);
      }
    }

    /**
     * Returns the mask representing the current estimation of the foreground.
     * @return Foregournd mask.
     */
    void getMask(ImageBool& imgBoolForeground) {
      if (MORPH) {
	imgBoolForeground = imgBoolMorph;
      } else {
	imgBoolForeground = imgBoolClean;
      }
    }

    /**
     * Returns the current background estimation.
     * @return Background estimation..
     */
    void getBackground(ImageGray8& imgBackground) {
      	imgBackground = imgSD_M;
    }

    private:
    ImageGray8 imgSD_M;
    ImageGray8 imgSD_D;
    ImageGray8 imgSD_V;
    ImageGray8 imgSD_U_M;
    ImageGray8 imgSD_U_STD;
    ImageBool  imgBoolMorph;
    ImageBool  imgBoolOut;
    ImageBool  imgBoolClean;
    Element elementEro;
    Element elementDil;

    
    /*********************************************************************
     * CleanImage()
     *
     * arguments:
     * output:
     *
     * Cleans the image by making the regions more clear.
     *********************************************************************/
    void CleanImage(ImageBool& imgIn, ImageBool& imgOut) {
      ImageBool::pixel_type pix1, pix2, pix_end;
      ImageBool::point_type posBool, offset;
      bool n, s, e, w;
  
      //Clean the image
      for(pix1=imgIn.begin(), pix_end=imgIn.end(), pix2=imgOut.begin();
	  pix1 != pix_end;
	  ++pix1, ++pix2) {
	if(*pix1) {
	  posBool = !pix2 + offset( 1, 0);
	  w = !(imgIn(posBool));
	  posBool = !pix2 + offset(-1, 0);
	  e = !(imgIn(posBool));
	  posBool = !pix2 + offset( 0, 1);
	  s = !(imgIn(posBool));
	  posBool = !pix2 + offset( 0,-1);
	  n = !(imgIn(posBool));
      
	  *pix2 = !((n&&s)||(e&&w));
	} else {
	  *pix2 = false;
	}
      }
    }



    /**
     * Removes shadows from the foreground estimation.
     * @param imgIn The next image in the series.
     * @paral imgOut The output mask with shadows removed.
     */
    void ShadowRemove(ImageGray8& imgIn,
		      ImageBool&  imgOut) {
      ImageBool::pixel_type         pix, pix_end;
      int                           j, k, temp, test, M;
      double                        dTemp, devTemp, devMean, devStd, Msize, Lstd, Llow;
      ImageGray8::point_type        pos;
      mirage::img::Coordinate       origin;
      
      M = PARAM::shadowWindow();
      
      Msize = (2*((double)M)+1) * (2*((double)M)+1);
      
      Lstd = PARAM::shadowVar();
      Llow = PARAM::shadowRatio();
      
      //Shadow detection loop
      for(pix = imgOut.begin(), pix_end = imgOut.end();
	  pix != pix_end; ++pix) {
	
	test = *pix;
	
	if (test > 0) {
	  //Foreground, do shadow detection
	  devTemp = 0;
	  devMean = 0;
	  devStd = 0;
	  origin = !pix;
	  for (j = (-1*M); j < (M+1); j++) {   //calculate std deviation
	    for (k = (-1*M); k < (M+1); k++) {
	      pos = origin + mirage::img::Coordinate(j, k);
	      dTemp = ((double) imgIn(pos)) / ((double) imgSD_M(pos));
	      devTemp = devTemp + dTemp*dTemp;
	      devMean = devMean + dTemp;
	    }
	  }
	  devMean = devMean / Msize;
	  devStd = sqrt(devTemp/Msize - devMean*devMean);
	  dTemp = ((double) imgIn(origin)) / ((double) imgSD_M(origin));
	  
	  if ((devStd < Lstd) && (dTemp < 1) && (dTemp >= Llow)) {
	    //Shadow!
	    *pix = 0;
	  }
	}
      }
    }

    //Copies the image, pixel-by-pixel
    class CopyImageOp {
    public:
    inline void operator()(mirage::colorspace::GRAY_8& src,
			   mirage::colorspace::GRAY_8& res) {
      res = (mirage::colorspace::GRAY_8) (src);
    }
    };
    typedef mirage::algo::UnaryOp<ImageGray8,ImageGray8,CopyImageOp> CopyImage;
    
    //Initialzes an image to a given value
    class InitImageOp {
    public:
      inline void operator()(mirage::colorspace::GRAY_8& res,
			     mirage::colorspace::GRAY_8& src) {
	res = (mirage::colorspace::GRAY_8) (PARAM::stdMin());
      }
    };
    typedef mirage::algo::UnaryOp<ImageGray8,ImageGray8,InitImageOp> InitImage;
    
    //Creates a blank black image
    class MakeNullOp {
    public:
      inline void operator()(mirage::colorspace::GRAY_8& res,
			     mirage::colorspace::GRAY_8& src) {
	res = (mirage::colorspace::GRAY_8) (0);
      }
    };
    typedef mirage::algo::UnaryOp<ImageGray8,ImageGray8,MakeNullOp> MakeNull;
    
    };

  }

}

#endif
