/*
    Nebulabrot implementation in C++
    adapted from Paul Nylander (http://nylander.wordpress.com, http://www.bugman123.com)
*/

#include <iostream>
#include <complex>
#include <vector>
#include <algorithm>
#include <dispatch/dispatch.h>
#include <gmp.h>
#include <mpc.h>
#include <mpfr.h>
#include <ext/malloc_allocator.h>

#include "pngwriter.h"

using namespace std;

#define FloatType long double
#define PRECISION 256
#define ROUNDMODE MPC_RNDDD 
#define FROUNDMODE GMP_RNDD 


//typedef complex<long double> ComplexDouble;
typedef mpc_t ComplexDouble;

struct Color {
    FloatType red;
    FloatType green;
    FloatType blue;
};

int Mandelbrot(ComplexDouble z) {
   
    mpc_t* value_list = (mpc_t*)malloc(1000);

    mpc_t y;
    mpc_init3(y, PRECISION, PRECISION);
    mpc_set(y, z, ROUNDMODE);

    mpc_init3(value_list[0], PRECISION, PRECISION);
    mpc_set_d_d(value_list[0], 0.0, 0.0, ROUNDMODE);

    mpfr_t absolute_value;
    mpfr_init2(absolute_value, 512);
    mpc_abs(absolute_value, value_list[0], FROUNDMODE);

    int i;
    for(i = 0; mpfr_get_ld(absolute_value, FROUNDMODE) < 2.0 && i<1000; i++) {
        //y = std::pow(list.back(), 2) + z;
	mpc_pow_ld(y, value_list[i], 2.0, ROUNDMODE);
	mpc_add(y, y, z, ROUNDMODE);
	mpc_init3(value_list[i+1], PRECISION, PRECISION);
	mpc_set(value_list[i+1], y, ROUNDMODE);
	mpc_abs(absolute_value, value_list[i+1], FROUNDMODE);
    }
    
    cout << i << endl;
    return i;
}

int CreateNebulabrotArray(mpc_t** image, int imageSize, FloatType xMin, FloatType xMax, FloatType yMin, FloatType yMax) {
    FloatType xIntervalStep = abs(xMin-xMax)/(double(imageSize));
    FloatType yIntervalStep = abs(yMin-yMax)/(double(imageSize));
    
    unsigned int iterations = 0, matches = 0;
    
    cout << "STATUS: Interval stepping ∆x, ∆y = " << xIntervalStep << ", " << yIntervalStep << endl;
   
#ifdef USE_GCD
    // TODO: change code in GMP-compatible manner
    dispatch_queue_t queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_DEFAULT, 0);
#endif

    for(FloatType x = xMin; x <= xMax; x += xIntervalStep) {	    
#ifdef USE_GCD
	// TODO: change code in GMP-compatible manner
	dispatch_apply((yMax-yMin)/yIntervalStep, queue, ^(size_t yPoint) {
	FloatType y = yMin + yPoint*yIntervalStep;
#else
	for(FloatType y = yMin; y <= yMax; y += yIntervalStep) {
#endif
	    mpc_t xy;
	    mpc_init3(xy, PRECISION, PRECISION);
	    mpc_set_ld_ld(xy, x, y, ROUNDMODE);
	    mpc_out_str(NULL, 10, 32, xy, ROUNDMODE);

            if(Mandelbrot(xy) < 1000) {
	    	mpc_t z;
		mpc_init3(z, PRECISION, PRECISION);
		mpc_set_ld_ld(z, 0.0, 0.0, ROUNDMODE);

		mpfr_t absolute_value;
		mpfr_init2(absolute_value, 2*PRECISION);
		mpc_abs(absolute_value, z, FROUNDMODE);
                while(mpfr_get_ld(absolute_value, FROUNDMODE) < 2.0) {
		    mpc_pow_ld(z, z, 2.0, ROUNDMODE);
		    mpc_add(z, z, xy, ROUNDMODE);
		                        
                    int i = floor((imageSize*mpfr_get_ld(mpc_imagref(z), FROUNDMODE)+1.5)/3);
                    int j = -floor((imageSize*mpfr_get_ld(mpc_realref(z), FROUNDMODE)+2)/3);
                    
                    i += floor(imageSize/2);
                    j -= floor(2*imageSize/3);
                    
                    if(i < 0) {
                        i = imageSize + i;
                    }
                    if(j < 0) {
                        j = imageSize + j;
                    }
                    
                    if(i < imageSize && j < imageSize && i > 0 && j > 0) {
			//mpc_add(image[i][j], image[i][j], xy, ROUNDMODE);
                        //matches++;
                    }
                }
            }
            //iterations++;
            //if(iterations%10000 == 0.0) {
            //    cout << iterations << ":" << matches << endl;
            //}
#ifdef USE_GCD
	});
#else
	}
#endif
    }
    
    return 0;
}

Color HueToRGB(FloatType h, FloatType s, FloatType l) {
    Color tVec;
    
    if(s == 0.0) {
        tVec.red = l;
        tVec.green = l;
        tVec.blue = l;
        return tVec;
    }
    
    FloatType q, p, hk, tRed, tGreen, tBlue;
    
    if(l < 0.5) {
        q = l*(1.0+s);
    } 
    if(l >= 0.5) {
        q = l+s-l*s;
    }
    
    hk = h/360.0;
    p = 2*l-q;
    
    tVec.red = hk + 1/3;
    tVec.green = hk;
    tVec.blue = hk - 1/3;
    
    FloatType* clr = &(tVec.red);
    for(int i=0; i<sizeof(Color); i+=sizeof(Color)/3) {
        clr += i;
        if(*clr < 0.0) {
            *clr = *clr + 1.0;
        }
        if(*clr > 1.0) {
            *clr = *clr - 1.0;
        }
        
        if((*clr)<1/6) {
            *clr = p + ((q-p)*6*(*clr));
        } if((*clr) >= 1/6 && (*clr) < 0.5) {
            *clr = q;
        } if((*clr) >= 0.5 && (*clr) < 2/3) {
            *clr = p + ((q-p)*6*(2/3-(*clr)));
        } else {
            *clr = p;
        }
    }
    
    return tVec;
}

void MapComplexToColor(mpc_t** source, Color** image, int imageSize) {
    pngwriter png(imageSize, imageSize, 0, "out.png");
    
    for(int i = 0; i<imageSize; i++) {
        for(int j = 0; j<imageSize; j++) {
            //Color clr = HueToRGB(std::arg(source[i][j])/(4.0*atan(1.0)), 1.0, min(std::abs(source[i][j])/18.0, 1.0));
	    mpfr_t argument;
	    mpfr_t absolute_value;
	    mpfr_init2(argument, 2*PRECISION);
	    mpfr_init2(absolute_value, 2*PRECISION);

	    mpc_arg(argument, source[i][j], FROUNDMODE);
	    mpc_abs(absolute_value, source[i][j], FROUNDMODE);

            png.plotHSV_blend(i, j, 1.0, 
			    mpfr_get_ld(argument, FROUNDMODE)/(1.0*(4.0*atan(1.0))), 
			    1.0, 
			    min(mpfr_get_ld(absolute_value, FROUNDMODE)/10.0, 1.0l)
			    );
            //png.plot(i, j, std::arg(source[i][j])/(4.0*atan(1.0)), min(std::abs(source[i][j])/18.0, 1.0), 1.0);
        }
    }
    
    png.setcompressionlevel(0); // To speed things up. Set to 0 for max speed when close()ing the image.
    png.setgamma(0.7);
    png.close();
}

int main(int argc, char** argv) {
    
    int size = int(strtol(argv[1], NULL, 10));

    cout << "STATUS: GMP " << __GNU_MP_VERSION << ", MPFR " << mpfr_get_version() << ", MPC " << mpc_get_version() << endl;
    
    cout << "STATUS: mpc_t size = " << sizeof(mpc_t) << endl;
    
    //mpc_t** image = (mpc_t**)malloc(size*sizeof(mpc_t*));
    //for(int i = 0; i < size; i++) {
    //    image[i] = (mpc_t*)malloc(size*sizeof(mpc_t));
    //}

    mpc_t** image = new mpc_t*[size];
    for(int i = 0; i < size; i++) {
        image[i] = new mpc_t[size];
    }


    cout << "STATUS: Now creating complex data array (" << size << ")..." << endl;
    CreateNebulabrotArray(image, size, -2.0, 1.0, -1.5, 1.5);
    
    Color** rgbData = new Color*[size];
    for(unsigned int i = 0; i < size; i++) {
          rgbData[i] = new Color[size];
    }
    
    cout << "STATUS: Now converting to RGB color space..." << endl;
    MapComplexToColor(image, rgbData, size);
        
    return 0;
}
