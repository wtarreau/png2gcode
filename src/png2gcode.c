/* Uses the simplified API, thus requires libpng 1.6 or above */
#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <png.h>
#include "font5x7.h"

/* default spindle level and feed rate in g-code mode. */
#define DEFAULT_SPINDLE_SHOW     1
#define DEFAULT_SPINDLE_RASTER   255
#define DEFAULT_FEED_SHOW        4800
#define DEFAULT_FEED_RASTER      1200

/* max number of arguments for a transformation */
#define MAX_XFRM_ARGS 10

/* 5x7 fonts give 6x8 as we keep one pixel around */
#define TEXT_HEIGHT 8
#define TEXT_WIDTH 6

struct image {
	uint32_t w;    // width in pixels
	uint32_t h;    // height in pixels
	uint8_t *rgba; // RGBA buffer
	uint8_t *gray; // grayscale buffer
	float *work;   // work area: 0.0 = white, 1.0+ = black
	float orgx;    // image horizontal origin in millimeters
	float orgy;    // image vertical origin in millimeters
	float mmw;     // image width in millimiters
	float mmh;     // image height in millimeters
	float diam;    // max computed diameter (if -cC) or zero
	uint8_t crop_threshold; // crop above this value
	uint32_t minx, maxx;
	uint32_t miny, maxy;
};

enum out_fmt {
	OUT_FMT_NONE = 0,
	OUT_FMT_PNG,
	OUT_FMT_GCODE,
};

enum out_center {
	OUT_CNT_NONE = 0,
	OUT_CNT_AXIS,
	OUT_CNT_CIRC,
};

enum opt {
	OPT_CROP_BOTTOM = 256,
	OPT_CROP_LEFT,
	OPT_CROP_RIGHT,
	OPT_CROP_TOP,
	OPT_AUTO_CROP,
	OPT_SOFTEN,
	OPT_IMGW,
	OPT_IMGH,
	OPT_IMGD,
	OPT_ORGX,
	OPT_ORGY,
	OPT_PIXW,
	OPT_PIXH,
	OPT_PIXS,
	OPT_RL_SHIFT,
	OPT_RL_DELAY,
	OPT_LASER_ON,
	OPT_TEST,
	OPT_MAP,
	OPT_MAX_LENGTH,
	OPT_EXT_BOTTOM,
	OPT_EXT_LEFT,
	OPT_EXT_RIGHT,
	OPT_EXT_TOP,
	OPT_TEXT_FG_ALPHA,
	OPT_TEXT_FG_VALUE,
	OPT_TEXT_BG_ALPHA,
	OPT_TEXT_BG_VALUE,
	OPT_TEXT_ZOOM,
	OPT_TEXT,
};

const struct option long_options[] = {
	{"help",        no_argument,       0, 'h'              },
	{"in",          required_argument, 0, 'i'              },
	{"out",         required_argument, 0, 'o'              },
	{"center",      required_argument, 0, 'c'              },
	{"fmt",         required_argument, 0, 'f'              },
	{"add",         required_argument, 0, 'a'              },
	{"offset",      required_argument, 0, 'O'              },
	{"mul",         required_argument, 0, 'm'              },
	{"gamma",       required_argument, 0, 'g'              },
	{"gamma2",      required_argument, 0, 'G'              },
	{"hash",        no_argument,       0, 'H'              },
	{"twins",       no_argument,       0, 't'              },
	{"vertical",    no_argument,       0, 'V'              },
	{"diffuse",     required_argument, 0, 'd'              },
	{"normalize",   no_argument,       0, 'n'              },
	{"raw-preview", no_argument,       0, 'r'              },
	{"quantize",    required_argument, 0, 'q'              },
	{"qfreq",       required_argument, 0, 'Q'              },
	{"mode",        required_argument, 0, 'M'              },
	{"feed",        required_argument, 0, 'F'              },
	{"passes",      required_argument, 0, 'P'              },
	{"spindle",     required_argument, 0, 'S'              },
	{"beam-width",  required_argument, 0, 'w'              },
	{"x-accel",     required_argument, 0, 'A'              },
	{"zero",        no_argument,       0, 'z'              },
	{"soften",      required_argument, 0, OPT_SOFTEN       },
	{"auto-crop",   required_argument, 0, OPT_AUTO_CROP    },
	{"crop-bottom", required_argument, 0, OPT_CROP_BOTTOM  },
	{"crop-left",   required_argument, 0, OPT_CROP_LEFT    },
	{"crop-right",  required_argument, 0, OPT_CROP_RIGHT   },
	{"crop-top",    required_argument, 0, OPT_CROP_TOP     },
	{"ext-bottom",  required_argument, 0, OPT_EXT_BOTTOM   },
	{"ext-left",    required_argument, 0, OPT_EXT_LEFT     },
	{"ext-right",   required_argument, 0, OPT_EXT_RIGHT    },
	{"ext-top",     required_argument, 0, OPT_EXT_TOP      },
	{"image-width", required_argument, 0, OPT_IMGW         },
	{"image-height",required_argument, 0, OPT_IMGH         },
	{"image-diam",  required_argument, 0, OPT_IMGD         },
	{"orig-x",      required_argument, 0, OPT_ORGX         },
	{"orig-y",      required_argument, 0, OPT_ORGY         },
	{"pixel-width", required_argument, 0, OPT_PIXW         },
	{"pixel-height",required_argument, 0, OPT_PIXH         },
	{"pixel-size",  required_argument, 0, OPT_PIXS         },
	{"rl-shift",    required_argument, 0, OPT_RL_SHIFT     },
	{"rl-delay",    required_argument, 0, OPT_RL_DELAY     },
	{"laser-on",    required_argument, 0, OPT_LASER_ON     },
	{"test",        required_argument, 0, OPT_TEST         },
	{"map",         required_argument, 0, OPT_MAP          },
	{"max-length",  required_argument, 0, OPT_MAX_LENGTH   },
	{"text-fg-alpha", required_argument, 0, OPT_TEXT_FG_ALPHA },
	{"text-fg-value", required_argument, 0, OPT_TEXT_FG_VALUE },
	{"text-bg-alpha", required_argument, 0, OPT_TEXT_BG_ALPHA },
	{"text-bg-value", required_argument, 0, OPT_TEXT_BG_VALUE },
	{"text-zoom",   required_argument, 0, OPT_TEXT_ZOOM    },
	{"text",        required_argument, 0, OPT_TEXT         },
	{0,             0,                 0, 0                }
};

/* describe one transformation to apply to the image */
enum xfrm_op {
	XFRM_NOP = 0,
	XFRM_ADD,       // args: 0=value
	XFRM_MUL,       // args: 0=value
	XFRM_GAM,       // args: 0=value
	XFRM_GAM2,      // args: 0=value
	XFRM_HASH,      // args: none
	XFRM_NORMALIZE, // args: none
	XFRM_QUANTIZE,  // args: 0=levels
	XFRM_SOFTEN,    // args: 0=value
	XFRM_QFREQ,     // args: 0=levels
	XFRM_TWINS,     // args: none
	XFRM_VERTICAL,  // args: none
	XFRM_MAP,       // args: [from_min:from_max:]to_min:to_max[:gamma2]
	XFRM_MAX_LENGTH,// args: 0=maxlen
	XFRM_OFST,      // args: 0=value
	XFRM_ZERO,      // args: none
	XFRM_TEXT,      // args: 0=xalign(-1,0,1) 1=yalign, 2=#lines, 3=zoom, 4=fgalpha, 5=fgcol, 6=bgalpha, 7=bgcol
};

struct xfrm {
	enum xfrm_op op;
	float args[MAX_XFRM_ARGS];
	char *text;
	struct xfrm *next;
};

/* type of operations that can be produced in a G-CODE pass */
enum pass_mode {
	PASS_MODE_ORIGIN = 0, // only move to origin
	PASS_MODE_X,          // move over the X axis
	PASS_MODE_Y,          // move over the Y axis
	PASS_MODE_AXIS,       // move over the X then the Y axis
	PASS_MODE_DIAG,       // move over the diagonal
	PASS_MODE_FRAME,      // move over the frame countour
	PASS_MODE_RASTER,     // raster image, bidirectional
	PASS_MODE_RASTER_LR,  // raster image left-to-right only
	PASS_MODE_CONTOUR,    // move along the image's countour
	PASS_MODES            // must be last one
};

const char *pass_mode_names[PASS_MODES] = {
	[PASS_MODE_ORIGIN] = "origin",
	[PASS_MODE_X]         = "x",
	[PASS_MODE_Y]         = "y",
	[PASS_MODE_AXIS]      = "axis",
	[PASS_MODE_DIAG]      = "diag",
	[PASS_MODE_FRAME]     = "frame",
	[PASS_MODE_RASTER]    = "raster",
	[PASS_MODE_RASTER_LR] = "raster-lr",
	[PASS_MODE_CONTOUR]   = "contour",
};

struct pass {
	enum pass_mode mode;
	int spindle;       // <0 = not set
	int feed;          // <0 = not set
	int passes;        // <=0 permitted (pass skipped)
	struct pass *next; // next pass or NULL
};

struct material {
	float diffusion;
} material = {
	      .diffusion = 0.18,  // 18% in each direction for each step results in 100 spread in 8 pixels
};

struct machine {
	float rl_shift;           // offset to add to R->L paths in raster mode to compensate for machine imprecision
	float rl_delay;           // time offset to add to R->L paths in raster mode to compensate for machine imprecision
	float beam_w;             // beam width in mm (limited to px size for now)
	float x_accel;            // acceleration/deceleration distance on X axis
	const char *laser_on;     // laser-on command
} machine = {
	     .rl_shift = 0,
	     .rl_delay = 0,
	     .beam_w   = 0,
	     .laser_on = "M4",
};

/* display the message and exit with the code */
__attribute__((noreturn)) void die(int code, const char *format, ...)
{
	va_list args;

	va_start(args, format);
	vfprintf(stderr, format, args);
	va_end(args);
	exit(code);
}

/* founds <f> to at most 4 digits after the point to reduce output size and
 * make sure %g will not emit "e-5" or so.
 */
double f4(double f)
{
	f = nearbyint(f * 10000.0) / 10000.0;
	return f;
}

/* returns the real X for an image from a dot position */
double imgxr(const struct image *img, int x)
{
	double xr;

	xr = x * img->mmw / img->w;
	xr = roundf(xr * 1000.0) / 1000.0;
	return xr;
}

/* returns the real Y for an image from a dot position */
double imgyr(const struct image *img, int y)
{
	double yr;

	yr = y * img->mmh / img->h;
	yr = roundf(yr * 1000.0) / 1000.0;
	return yr;
}

/* returns the pixel intensity at <x,y> or 0 if outside of the viewing area */
float get_pix(const struct image *img, uint32_t x, uint32_t y)
{
	if (y >= img->h)
		return 0;
	if (x >= img->w)
		return 0;
	return img->work[y * img->w + x];
}

/* computes the square of the distance between points 1 and 2 */
float sqdist(uint32_t x1, uint32_t y1, uint32_t x2, uint32_t y2, float resx, float resy)
{
	resx *= (int)(x2 - x1);
	resy *= (int)(y2 - y1);
	return resx * resx + resy * resy;
}

void usage(int code, const char *cmd)
{
	die(code,
	    "\n"
	    "Usage: %s [options]* [transformations]* [passes]*\n"
	    "  -h --help                    show this help\n"
	    "  -f --fmt <format>            output format (png, gcode), defaults to file ext\n"
	    "  -r --raw-preview             make the PNG output match exactly the laser dot intensity\n"
	    "  -i --in <file>               input PNG file name\n"
	    "     --auto-crop <threshold%%>  automatic crop of input image (suggested: 95%%)\n"
	    "     --crop-bottom <size>      crop this number of pixels from the bottom\n"
	    "     --crop-left   <size>      crop this number of pixels from the left\n"
	    "     --crop-right  <size>      crop this number of pixels from the right\n"
	    "     --crop-top    <size>      crop this number of pixels from the top\n"
	    "     --ext-bottom  <size>      extend by this number of pixels on the bottom\n"
	    "     --ext-left    <size>      extend by this number of pixels on the left\n"
	    "     --ext-right   <size>      extend by this number of pixels on the right\n"
	    "     --ext-top     <size>      extend by this number of pixels on the top\n"
	    "     --test <size>             generate a test pattern this size instead of input\n"
	    "  -o --out <file>              output file name\n"
	    "     --image-width  <size>     image width in millimeters (default: automatic)\n"
	    "     --image-height <size>     image height in millimeters (default: automatic)\n"
	    "     --image-diam   <size>     image diameter for Circle mode (default: automatic)\n"
	    "     --orig-x       <pos>      X origin position in millimeters (default: 0.0)\n"
	    "     --orig-y       <pos>      Y origin position in millimeters (default: 0.0)\n"
	    "     --pixel-width  <size>     pixel width in millimeters\n"
	    "     --pixel-height <size>     pixel height in millimeters\n"
	    "     --pixel-size   <size>     pixel size in millimeters (sets width and height)\n"
	    "  -c, --center <mode>          auto-center output coordinates (N=none, A=axis, C=circle)\n"
	    "Transformations are series of operations applied to the work area:\n"
	    "  -a --add <value>             add <value> [-1..1] to the intensity\n"
	    "  -O --offset <value>          add <value> [-1..1] to a non-zero intensity\n"
	    "  -z --zero                    turn lowest value to zero\n"
	    "  -g --gamma <value>           apply abs value from bottom if >0 or from top if <0 (def 1.0)\n"
	    "  -G --gamma2 <value>          apply a centered gamma correction (def 1.0)\n"
	    "  -H --hash                    hash the image by setting 50%% of the dots to intensity 0\n"
	    "  -t --twins                   average adjacent pixels and send as opposed twins\n"
	    "  -V --vertical                vertical stripes (optimal for cardboard)\n"
	    "  -m --mul <value>             multiply intensity by <value>\n"
	    "  -n --normalize               normalize work from 0.0 to 1.0\n"
	    "  -q --quantize <levels>       quantize to <levels> levels\n"
	    "  -Q --qfreq <levels>          quantize to <levels> levels based on frequency\n"
	    "     --soften <value>          subtract neighbors' average times <value>\n"
	    "     --map <args[2..5]>        map values [from_min:from_max:]to_min:to_max[:gamma2]\n"
	    "     --max-length <pixels>     limit long series for use with FS dithering (def: 0)\n"
	    "     --text-fg-value <value>   set FG value for future text (def: 1.0)\n"
	    "     --text-fg-alpha <value>   set FG alpha for future text (def: 1.0)\n"
	    "     --text-bg-value <value>   set BG value for future text (def: 0.0)\n"
	    "     --text-bg-alpha <value>   set BG alpha for future text (def: 0.0)\n"
	    "     --text-zoom <int>         set integer zoom value for future text (def: 1)\n"
	    "     --text [@<align>:]<msg>   write/append <msg>. align=[{bct}][{lcr}]\n"
	    "Passes are used in G-CODE output format (-f gcode):\n"
	    "  -M --mode    <mode>          pass mode (origin,x,y,axis,diag,frame,raster,raster-lr)\n"
	    "  -F --feed    <value>         feed rate (mm/min, def:4800 frame, 1200 raster)\n"
	    "  -S --spindle <value>         spindle value for intensity 1.0 (def:1 or 255)\n"
	    "  -P --passes  <value>         number of passes with previous parameters (def:1)\n"
	    "Material characteristics:\n"
	    "  -d --diffuse <ratio>         adjacent radiation for hash & soften (def:0.18)\n"
	    "Machine settings:\n"
	    "  -w --beam-width <mm>         beam width in millimeters (def: 0, max:90%% pxw, neg=alt)\n"
	    "  -A --x-accel <mm>            X axis acceleration distance in millimeters (def: 0)\n"
	    "     --rl-shift <millimeters>  offset to apply to R->L path in raster mode (def:0)\n"
	    "     --rl-delay <millisec>     time offset to apply to R->L path in raster mode (def:0)\n"
	    "     --laser-on <cmd>          command to turn laser ON (def:M4)\n"
	    "Notes:\n"
	    "  - for images, use -H for wood or -t on aluminum\n"
	    "  - an empty image may be created from just text\n"
	    "  - for pure text, -d0 is recommended to sharpen the text\n"
	    "\n", cmd);
}

/* reads file <file> into rgba image <img>, which will be initialized. The image
 * will go from bottom to top to accommodate from GCODE's image directions, but
 * this can be changed by setting the row_stride argument to 1 instead of -1.
 * Returns non-zero on sucess, or 0 on failure, in which case the contents of
 * <img> remains undefined.
 */
int read_rgba_file(const char *file, struct image *img)
{
	const int row_stride = -1; // bottom to top
	png_image png_image;

	memset(img, 0, sizeof(*img));
	memset(&png_image, 0, sizeof(png_image));

	png_image.version = PNG_IMAGE_VERSION;
	if (!png_image_begin_read_from_file(&png_image, file))
		goto fail;

	/* process all images as RGBA so that we can turn transparent into white*/
	png_image.format  = PNG_FORMAT_RGBA;

	img->crop_threshold = 255;
	img->w = png_image.width;
	img->h = png_image.height;
	img->rgba = malloc(PNG_IMAGE_SIZE(png_image));
	if (!img->rgba)
		goto fail;

	if (!png_image_finish_read(&png_image, NULL, img->rgba, row_stride * PNG_IMAGE_ROW_STRIDE(png_image), NULL))
		goto fail;

	return 1;

 fail:
	free(img->rgba);
	return 0;
}

/* generate a test pattern <w> pixels wide and <h> pixels high.
 * Returns non-zero on sucess, or 0 on failure, in which case the contents of
 * <img> remains undefined.
 */
int generate_test_pattern(struct image *img, uint32_t w, uint32_t h)
{
	uint32_t x, y, v;

	memset(img, 0, sizeof(*img));

	img->crop_threshold = 255;
	img->w = w;
	img->h = h;

	img->gray = malloc(img->w * img->h * sizeof(*img->gray));
	if (!img->gray)
		return 0;

	for (y = 0; y < img->h; y++) {
		for (x = 0; x < img->w; x++) {
			v = (x ^ y) & 1;
			v *= 255;
			img->gray[y * img->w + x] = v;
		}
	}

	img->minx = 0;
	img->maxx = img->w - 1;
	img->miny = 0;
	img->maxy = img->h - 1;

	return 1;
}

/* free all buffers in image <img> */
void free_image(struct image *img)
{
	free(img->rgba);
	free(img->gray);
	free(img->work);
	memset(img, 0, sizeof(*img));
}

/* write the gray buffer from <img> into file <file>, or to stdout if <file> is
 * NULL. The image will go from bottom to top to accommodate from GCODE's image
 * directions, but this can be changed by setting the row_stride argument to 1
 * instead of -1. Returns non-zero on success, otherwise zero.
 */
int write_gray_file(const char *file, struct image *img)
{
	const int row_stride = -1; // bottom to top
	png_image png_image;
	int ret;

	memset(&png_image, 0, sizeof(png_image));
	png_image.version = PNG_IMAGE_VERSION;
	png_image.width   = img->w;
	png_image.height  = img->h;
	png_image.format  = PNG_FORMAT_GRAY;

	if (file)
		ret = png_image_write_to_file(&png_image, file, 0, img->gray, row_stride * PNG_IMAGE_ROW_STRIDE(png_image), NULL);
	else
		ret = png_image_write_to_stdio(&png_image, stdout, 0, img->gray, row_stride * PNG_IMAGE_ROW_STRIDE(png_image), NULL);
	return ret;
}

/* convert image <img> from rgba to gray. Transparent is turned to white. The
 * gray buffer will be allocated and the rgba buffer will be freed. Non-zero is
 * returned on success, zero on error.
 */
int rgba_to_gray(struct image *img)
{
	uint32_t x, y;
	uint32_t minx, maxx, miny, maxy;

	if (!img->rgba)
		return 0;

	img->gray = malloc(img->w * img->h * sizeof(*img->gray));
	if (!img->gray)
		return 0;

	maxx = 0; minx = img->w - 1;
	maxy = 0; miny = img->h - 1;

	for (y = 0; y < img->h; y++) {
		for (x = 0; x < img->w; x++) {
			uint32_t v, p;
			uint8_t r, g, b, a;

			p = (y * img->w + x) * 4;
			r = img->rgba[p + 0];
			g = img->rgba[p + 1];
			b = img->rgba[p + 2];
			a = img->rgba[p + 3];
			v = a * (r + g + b) + (255 - a) * (3 * 255);
			v /= 3 * 255;
			img->gray[y * img->w + x] = v;
			if (v <= img->crop_threshold) {
				if (x < minx)
					minx = x;
				if (x > maxx)
					maxx = x;
				if (y < miny)
					miny = y;
				if (y > maxy)
					maxy = y;
			}
		}
	}

	/* if no pixel was found, do not crop */
	if (minx > maxx) {
		minx = 0;
		maxx = img->w - 1;
	}

	if (miny > maxy) {
		miny = 0;
		maxy = img->h - 1;
	}

	img->minx = minx;
	img->maxx = maxx;
	img->miny = miny;
	img->maxy = maxy;

	free(img->rgba);
	img->rgba = NULL;
	return 1;
}

/* convert image <img> from gray to work. White is turned to strength 0.0, and
 * black is turned to strength 1.0. The work buffer will be allocated and the
 * gray buffer will be freed. Non-zero is returned on success, zero on error.
 */
int gray_to_work(struct image *img)
{
	uint32_t x, y;

	if (!img->gray)
		return 0;

	img->work = malloc(img->w * img->h * sizeof(*img->work));
	if (!img->work)
		return 0;

	for (y = 0; y < img->h; y++) {
		for (x = 0; x < img->w; x++) {
			uint32_t p = y * img->w + x;

			img->work[p] = (255 - img->gray[p]) / 255.0;
		}
	}

	free(img->gray);
	img->gray = NULL;
	return 1;
}

/* convert image <img> from work to gray. Signal strengths 0.0 and lower are
 * turned to white, and strengths 1.0 and higher are turned to black. The gray
 * buffer will be allocated and the work buffer will be freed. Non-zero is
 * returned on success, zero on error. If <raw> is non-zero, the exact laser
 * dot intensity is reported for each pixel, otherwise the previous considers
 * the material's diffusion and tries to provide an estimate of what the output
 * would look like.
 */
int work_to_gray(struct image *img, int raw)
{
	uint32_t x, y;

	if (!img->work)
		return 0;

	img->gray = malloc(img->w * img->h * sizeof(*img->gray));
	if (!img->gray)
		return 0;

	for (y = 0; y < img->h; y++) {
		for (x = 0; x < img->w; x++) {
			uint32_t p = y * img->w + x;
			float    v = img->work[p];
			float d = material.diffusion;
			float d2 = 2*d*d;

			/* quadratic distance effect hence 2*d^2 for diagonal
			 * since it spreads through two adjacent pixels. Also
			 * don't count current pixel.
			 */
			if (!raw) {
				v += get_pix(img, x - 1, y - 1) * d2;
				v += get_pix(img, x + 0, y - 1) * d;
				v += get_pix(img, x + 1, y - 1) * d2;

				v += get_pix(img, x - 1, y) * d;
				v += get_pix(img, x + 1, y) * d;

				v += get_pix(img, x - 1, y + 1) * d2;
				v += get_pix(img, x + 0, y + 1) * d;
				v += get_pix(img, x + 1, y + 1) * d2;
			}
			v = 255 * (1.0 - v);
			if (v < 0)
				v = 0;
			else if (v > 255.0)
				v = 255.0;
			img->gray[p] = (uint8_t)v;
		}
	}

	free(img->work);
	img->work = NULL;
	return 1;
}

/* crop a gray image <img>, to keep only (<x0>,<y0)-(<x1>,<y1>), all included.
 * First columns and rows are numbered zero. Returns non-zero on success, 0 on
 * failure. x0, x1, y0, y1 must be within the original buffer dimensions, with
 * x0<=x1, y0<=y1. The image's w and h are updated. The buffer is rearranged
 * and the extra area is released. Note that the img->area pointer may change.
 */
int crop_gray_image(struct image *img, uint32_t x0, uint32_t y0, uint32_t x1, uint32_t y1)
{
	int row_pre, row_post;
	uint8_t *src, *dst;
	uint32_t x, y;

	if (x0 >= img->w || x1 >= img->w || x0 > x1)
		return 0;

	if (y0 >= img->h || y1 >= img->h || y0 > y1)
		return 0;

	if (!img->gray)
		return 0;

	row_pre = x0;
	row_post = img->w - 1 - x1;

	src = dst = img->gray;
	src += y0 * img->w;
	for (y = y0; y <= y1; y++) {
		src += row_pre;
		for (x = x0; x <= x1; x++)
			*dst++ = *src++;
		src += row_post;
	}

	img->w = x1 - x0 + 1;
	img->h = y1 - y0 + 1;

	dst = realloc(img->gray, img->w * img->h * sizeof(*img->gray));
	if (!dst)
		return 0;

	img->gray = dst;
	return 1;
}

/* extend gray image <img>, by <l> pixels on the left, <r> pixels on the right,
 * <t> pixels at the top, <b> pixels at the bottom. The area will be
 * reallocated if needed. Returns non-zero on success, 0 on failure. The
 * image's w and h are updated.
 */
int extend_gray_image(struct image *img, int l, int r, int t, int b)
{
	uint8_t *src, *dst;
	int x, y;
	int nw, nh;
	int fill = 255; // fill color = white

	nw = l + img->w + r;
	nh = t + img->h + b;
	dst = realloc(img->gray, nw * nh * sizeof(*img->gray));
	if (!dst)
		return 0;

	img->gray = dst;

	/* The image remains packed at the top, so let's start from the end
	 * to add the missing pixels.
	 */
	if (t)
		memset(img->gray + (b + img->h) * nw * sizeof(*img->gray), fill, t * nw * sizeof(*img->gray));

	/* process existing lines, possibly padded */
	for (y = b + img->h - 1; y >= b; y--) {
		src = img->gray + (y - b) * img->w + (img->w - 1);
		dst = img->gray + y * nw + nw - 1;
		// right padding
		for (x = nw - 1; x >= (int)(l + img->w); x--)
			*dst-- = fill;
		// copy
		for (; x >= l; x--)
			*dst-- = *src--;
		// left padding
		for (; x >= 0; x--)
			*dst-- = fill;
	}

	/* optionally pad the top*/
	if (b)
		memset(img->gray, fill, b * nw * sizeof(*img->gray));

	img->w = nw;
	img->h = nh;
	return 1;
}

/* append a transformation to list <curr> */
struct xfrm *xfrm_new(struct xfrm *curr, enum xfrm_op op, unsigned nbargs, const float *args)
{
	struct xfrm *new;

	if (nbargs > MAX_XFRM_ARGS)
		return NULL;

	new = calloc(1, sizeof(*new));
	if (!new)
		return NULL;

	new->op    = op;
	if (nbargs)
		memcpy(new->args, args, nbargs * sizeof(*args));

	if (curr)
		curr->next = new;
	return new;
}

/* append a text-based transformation to list <curr> */
struct xfrm *xfrm_new_text(struct xfrm *curr, enum xfrm_op op, const char *msg, unsigned zoom, float fg_alpha, float fg_value, float bg_alpha, float bg_value)
{
	int xalign = 0, yalign = 0;
	struct xfrm *new = curr;
	const char *p = msg;
	int lines = 1;

	if (curr && curr->op == op) {
		xalign = curr->args[0];
		yalign = curr->args[1];
	}

	if (*p == '@') {
		switch (p[1]) {
		case 'b': yalign = -1; p++; break;
		case 'c': yalign =  0; p++; break;
		case 't': yalign =  1; p++; break;
		}
		switch (p[1]) {
		case 'l': xalign = -1; p++; break;
		case 'c': xalign =  0; p++; break;
		case 'r': xalign =  1; p++; break;
		}
		while (p[1] != ':')
			p++;
		if (p[1] == ':')
			p++;
		msg = p + 1;
	}

	if (!curr || op != curr->op ||
	    xalign != curr->args[0] || yalign != curr->args[1] || zoom != curr->args[3] ||
	    fg_alpha != curr->args[4] || fg_value != curr->args[5] ||
	    bg_alpha != curr->args[6] || bg_value != curr->args[7]) {
		new = calloc(1, sizeof(*new));
		if (!new)
			return NULL;

		if (curr)
			curr->next = new;
	}

	/* count lines */
	p = msg - 1;
	while ((p = strchr(p + 1, '\n')))
		lines++;

	/* fill/update xfrm */
	new->op       = op;
	new->args[0]  = xalign;
	new->args[1]  = yalign;
	new->args[2] += lines;
	new->args[3]  = zoom;
	new->args[4]  = fg_alpha;
	new->args[5]  = fg_value;
	new->args[6]  = bg_alpha;
	new->args[7]  = bg_value;
	if (!new->text) {
		new->text = strdup(msg);
	}
	else {
		size_t size = strlen(new->text) + 1 + strlen(msg) + 1;
		char *txt = malloc(size);

		if (!txt) {
			if (new != curr)
				free(new);
			return NULL;
		}
		snprintf(txt, size, "%s\n%s", new->text, msg);
		free(new->text);
		new->text = txt;
	}
	return new;
}

/* Apply text described in <xfrm> to image <img> */
void xfrm_text(struct image *img, struct xfrm *xfrm)
{
	char *line, *next;
	int lines, linenum;
	int xalign, yalign, zoom;
	float fg_alpha, fg_value, bg_alpha, bg_value;
	int xinit, yinit, xchar, x, y;
	int width;

	if (xfrm->op != XFRM_TEXT)
		return;

	xalign = xfrm->args[0];
	yalign = xfrm->args[1];
	lines = xfrm->args[2];
	zoom = xfrm->args[3];
	fg_alpha = xfrm->args[4];
	fg_value = xfrm->args[5];
	bg_alpha = xfrm->args[6];
	bg_value = xfrm->args[7];

	if (yalign < 0) // bottom
		yinit = lines * zoom * TEXT_HEIGHT - 1;
	else if (yalign == 0) // center
		yinit = img->h / 2 + lines * zoom * TEXT_HEIGHT / 2 - 1;
	else // top
		yinit = img->h - 1;

	if (xalign < 0) // left
		xinit = 0;
	else if (xalign == 0) // center (will be adjusted for each line)
		xinit = img->w / 2;
	else // right
		xinit = img->w; // text length will be subtracted

	line = next = xfrm->text;
	for (linenum = 0; linenum < lines; linenum++) {
		next = strchr(line, '\n');
		if (!next)
			next = line + strlen(line);
		if (yinit < zoom * TEXT_HEIGHT - 1 || yinit >= img->h)
			goto skip;

		/* we're going to print between line and next-1 */
		width = next - line;
		if (xalign < 0)
			xchar = 0;
		else if (xalign == 0)
			xchar = xinit - width * zoom * TEXT_WIDTH / 2;
		else
			xchar = xinit - width * zoom * TEXT_WIDTH;

		while (xchar < 0) {
			xchar += zoom * TEXT_WIDTH;
			line++;
		}

		while (line < next) {
			char c = *line++;

			if (xchar > img->w - zoom * TEXT_WIDTH)
				break;

			if (c < 0x20 || c > 0x7f)
				c = 0x20;
			c -= 0x20;

			for (y = 0; y < zoom * TEXT_HEIGHT; y++) {
				if (y >= zoom * (TEXT_HEIGHT - 1) && !linenum)
					break;
				for (x = 0; x < zoom * TEXT_WIDTH; x++) {
					if (x >= zoom * (TEXT_WIDTH - 1) && line >= next)
						break;
					if (x < zoom * (TEXT_WIDTH - 1) && y < zoom * (TEXT_HEIGHT - 1) &&
					    font5x7[(uint8_t)c][y / zoom] & (1 << ((TEXT_WIDTH - 2) - (x / zoom))))
						img->work[(yinit - y) * img->w + xchar + x] =
							img->work[(yinit - y) * img->w + xchar + x] * (1.0 - fg_alpha) + fg_alpha * fg_value;
					else
						img->work[(yinit - y) * img->w + xchar + x] =
							img->work[(yinit - y) * img->w + xchar + x] * (1.0 - bg_alpha) + bg_alpha * bg_value;
				}
			}
			xchar += zoom * TEXT_WIDTH;
		}
	skip:
		yinit -= zoom * TEXT_HEIGHT;
		if (*next)
			next++;
		line = next;
	}
}

/* try to figure the center of the smallest circle including all pixels,
 * returns it as well as the computed radius in non-null pointers. Returns
 * non-zero on success.
 */
int find_center(struct image *img, float *outx, float *outy, float *outrad)
{
	int *minx, *maxx, *miny, *maxy;
	int x, y, x2, y2, pt, npt;
	uint64_t sumx, sumy;
	int ctrx, ctry;
	int moves;
	int ret;
	float *radius;
	float rad, maxrad;
	float resx, resy;
	struct point {
		int x, y;
	} *pts;

	ret = 0;
	minx = malloc(img->h * sizeof(*minx));
	maxx = malloc(img->h * sizeof(*maxx));
	miny = malloc(img->w * sizeof(*miny));
	maxy = malloc(img->w * sizeof(*maxy));
	pts = calloc((img->w + img->h) * 2, sizeof(*pts));
	radius = calloc(img->w * img->h, sizeof(*pts));

	if (!minx || !maxx || !miny || !maxy || !pts || !radius)
		goto done;

	resx = img->mmh / img->h;
	resy = img->mmw / img->w;

	if ((!resx || !resy) && img->diam)
		resx = resy = 1.0;

	for (y = 0; y < img->h; y++) {
		minx[y] = maxx[y] = -1;
		for (x = 0; x < img->w; x++)
			if (img->work[y * img->w + x] > 0.0) {
				minx[y] = x;
				break;
			}
		for (x = img->w - 1; x >= 0; x--)
			if (img->work[y * img->w + x] > 0.0) {
				maxx[y] = x;
				break;
			}
	}

	for (x = 0; x < img->w; x++) {
		miny[x] = maxy[x] = -1;
		for (y = 0; y < img->h; y++)
			if (img->work[y * img->w + x] > 0.0) {
				miny[x] = y;
				break;
			}
		for (y = img->h - 1; y >= 0; y--)
			if (img->work[y * img->w + x] > 0.0) {
				maxy[x] = y;
				break;
			}
	}

	/* store all boundaries */
	npt = 0;
	for (y = 0; y < img->h; y++) {
		if (minx[y] == -1)
			continue;
		pts[npt].x = minx[y];
		pts[npt].y = y;
		npt++;
		if (minx[y] == maxx[y])
			continue;
		pts[npt].x = maxx[y];
		pts[npt].y = y;
		npt++;
	}
	if (!npt)
		goto done;

	/* find an initial center. With this choice we know there are at
	 * least pixels in each direction.
	 */
	sumx = sumy = 0;
	for (pt = 0; pt < npt; pt++) {
		sumx += pts[pt].x;
		sumy += pts[pt].y;
	}
	ctrx = sumx / npt;
	ctry = sumy / npt;

	/* oscillate around the current center and move it towards the smallest max */
	moves = 0;
	while (1) {
		maxrad = 0;
		for (y = ctry - 1; y <= ctry + 1; y++) {
			if (y < 0 || y >= img->h)
				continue;

			for (x = ctrx - 1; x <= ctrx + 1; x++) {
				if (x < 0 || x >= img->w)
					continue;

				if (radius[y * img->w + x] > 0)
					continue; // already calculated
				maxrad = 0;
				for (pt = 0; pt < npt; pt++) {
					rad = sqdist(x, y, pts[pt].x, pts[pt].y, resx, resy);
					if (rad > maxrad)
						maxrad = rad;
				}
				radius[y * img->w + x] = maxrad;
			}
		}

		/* we have 9 measures around the current center, let's move in
		 * the most optimal direction now to pick the smallest radius,
		 * or stop if already at the center.
		 */

		x2 = ctrx; y2 = ctry;
		maxrad = 0;
		for (y = ctry - 1; y <= ctry + 1; y++) {
			if (y < 0 || y >= img->h)
				continue;

			for (x = ctrx - 1; x <= ctrx + 1; x++) {
				if (x < 0 || x >= img->w)
					continue;
				if (maxrad == 0 || radius[y * img->w + x] < maxrad) {
					maxrad = radius[y * img->w + x];
					x2 = x;
					y2 = y;
				}
			}
		}

		if (x2 == ctrx && y2 == ctry)
			break;

		moves++;
		ctrx = x2; ctry = y2;
	}

	rad = sqrt(radius[ctry * img->w + ctrx]);

	/* diameter was forced */
	if (img->diam) {
		img->mmw = img->w * img->diam / (rad * 2);
		img->mmh = img->h * img->diam / (rad * 2);
		resx *= img->diam / (rad * 2);
		resy *= img->diam / (rad * 2);
		rad   = img->diam / 2;
	}

	printf("Found optimal center at %d,%d in %d moves, max radius=%.2f mm, diam=%.2f mm\n",
	       ctrx, ctry, moves, rad, 2*rad);

	if (outx)
		*outx = ctrx * resx;

	if (outy)
		*outy = ctry * resy;

	if (outrad)
		*outrad = rad;
	ret = 1;

#ifdef DEBUG_CENTER
	/* hide all pixels that are not boundaries */
	for (y = 0; y < img->h; y++) {
		if (minx[y] == -1)
			continue;
		for (x = minx[y] + 1; x <= maxx[y] - 1; x++)
			if (y > miny[x] && y < maxy[x])
				img->work[y * img->w + x] = 0;
	}

	/* mark remaining pixels (2 per x, 2 per y) */
	for (y = 0; y < img->h; y++) {
		for (x = 0; x < img->w; x++)
			if (img->work[y * img->w + x] > 0) {
				img->work[y * img->w + x] = 1.0;
			}
	}

	/* mark the initial center */
	img->work[ctry * img->w + ctrx] = 0.5;
#endif

 done:
	free(radius);
	free(pts);
	free(maxy); free(miny);
	free(maxx); free(minx);
	return ret;
}

/* Measure the max width and height of text by cumulating L/C/R and T/C/B, by
 * scanning all xfrm_text.
 */
void measure_text_size(struct xfrm *xfrm, int *width, int *height)
{
	int l_w, c_w, r_w;
	int t_h, c_h, b_h;
	int w, wmax, h;
	char *line, *next;

	*width = *height = 0;
	l_w = c_w = r_w = t_h = c_h = b_h = 0;
	while (xfrm) {
		if (xfrm->op == XFRM_TEXT) {
			h = xfrm->args[2] * xfrm->args[3] * TEXT_HEIGHT; // lines * zoom
			if (xfrm->args[1] < 0 && b_h < h)       // yalign == bottom
				b_h = h;
			else if (xfrm->args[1] == 0 && c_h < h) // yalign == center
				c_h = h;
			else if (xfrm->args[1] > 0 && t_h < h)  // yalign == top
				t_h = h;

			line = next = xfrm->text;
			for (w = wmax = 0; *line; line = next) {
				next = strchr(line, '\n');
				if (!next)
					next = line + strlen(line);
				w = next - line;
				if (w > wmax)
					wmax = w;
				if (*next)
					next++;
				line = next;
			}

			w = wmax * xfrm->args[3] * TEXT_HEIGHT; // cols * zoom
			if (xfrm->args[0] < 0 && l_w < w)       // xalign == left
				l_w = w;
			else if (xfrm->args[0] == 0 && c_w < w) // xalign == center
				c_w = w;
			else if (xfrm->args[0] > 0 && r_w < w)  // xalign == right
				r_w = w;
		}
		xfrm = xfrm->next;
	}

	/* OK let's write this back now */
	*width = l_w + c_w + r_w;
	*height = b_h + c_h + t_h;
}
/* apply transformations starting at <xfrm> to image <img>. Return non-zero on
 * success, zero on failure.
 */
int xfrm_apply(struct image *img, struct xfrm *xfrm)
{
	uint32_t x, y;
	uint32_t *freq = NULL;
	float *soften = NULL;
	uint32_t q;
	float bottom = 0.0;

	while (xfrm) {
		float norm_min = 1.0;
		float norm_max = 0.0;

		if (xfrm->op == XFRM_SOFTEN || xfrm->op == XFRM_HASH) {
			soften = malloc(img->h * img->w * sizeof(*soften));
			if (!soften)
				return 0;

			for (y = 0; y < img->h; y++) {
				for (x = 0; x < img->w; x++) {
					float v = 0.0;
					float d = material.diffusion;
					float d2 = 2*d*d;

					/* quadratic distance effect hence 2*d^2 for diagonal
					 * since it spreads through two adjacent pixels. Also
					 * don't count current pixel.
					 */
					v += get_pix(img, x - 1, y - 1) * d2;
					v += get_pix(img, x + 0, y - 1) * d;
					v += get_pix(img, x + 1, y - 1) * d2;

					v += get_pix(img, x - 1, y) * d;
					v += get_pix(img, x + 1, y) * d;

					v += get_pix(img, x - 1, y + 1) * d2;
					v += get_pix(img, x + 0, y + 1) * d;
					v += get_pix(img, x + 1, y + 1) * d2;

					soften[y * img->w + x] = v;
				}
			}
		}
		else if (xfrm->op == XFRM_TWINS || xfrm->op == XFRM_VERTICAL) {
			/* we need the average of two adjacent X pixels */
			soften = malloc(img->h * img->w * sizeof(*soften));
			if (!soften)
				return 0;

			for (y = 0; y < img->h; y++)
				for (x = 0; x < img->w; x++) {
					float v = 0.0;
					float d = material.diffusion;
					float d2 = 2*d*d;

					/* quadratic distance effect hence 2*d^2 for diagonal
					 * since it spreads through two adjacent pixels. Also
					 * don't count current pixel.
					 */
					v += get_pix(img, x - 1, y - 1) * d2;
					v += get_pix(img, x + 0, y - 1) * d;
					v += get_pix(img, x + 1, y - 1) * d2;

					v += get_pix(img, x - 1, y) * d;
					v += get_pix(img, x + 1, y) * d;

					v += get_pix(img, x - 1, y + 1) * d2;
					v += get_pix(img, x + 0, y + 1) * d;
					v += get_pix(img, x + 1, y + 1) * d2;

					soften[y * img->w + x] = get_pix(img, x & ~1, y) + get_pix(img, x | 1, y) - v / 2;
				}
		}
		else if (xfrm->op == XFRM_NORMALIZE) {
			for (y = 0; y < img->h; y++) {
				for (x = 0; x < img->w; x++) {
					float v = img->work[y * img->w + x];
					if (v < norm_min)
						norm_min = v;
					if (v > norm_max)
						norm_max = v;
				}
			}

		}
		else if (xfrm->op == XFRM_QFREQ) {
			/* normalize values on 10 bits and store their frequencies */
			freq = calloc(1024, sizeof(*freq));
			if (!freq)
				return 0;

			for (y = 0; y < img->h; y++) {
				for (x = 0; x < img->w; x++) {
					uint32_t p = y * img->w + x;
					float    v = img->work[p];

					if (v <= 0.0)
						q = 0;
					else if (v >= 1.0)
						q = 1023;
					else
						q = 1023 * v;
					freq[q]++;
				}
			}

			/* make each value reflect its position in the spectrum.
			 * The last value will have the total number of pixels
			 * and should equal h*w. As such, now any pixel value v
			 * passed through freq[] can have its frequency normalized
			 * between [0..1] by dividing by w*h.
			 */
			for (q = 1; q < 1024; q++) {
				freq[q] += freq[q - 1];
			}
			//for (q = 0; q < 1024; q++) {
			//	printf("%u: %u (%f)\n", q, freq[q], (double)freq[q]/freq[1023]);
			//}
		}
		else if (xfrm->op == XFRM_ZERO) {
			/* find the lowest value */
			bottom = 1.0;

			for (y = 0; y < img->h; y++) {
				for (x = 0; x < img->w; x++) {
					uint32_t p = y * img->w + x;
					float    v = img->work[p];
					if (v < bottom)
						bottom = v;
				}
			}
		}
		else if (xfrm->op == XFRM_TEXT) {
			xfrm_text(img, xfrm);
		}

		for (y = 0; y < img->h; y++) {
			for (x = 0; x < img->w; x++) {
				uint32_t p = y * img->w + x;
				float    v = img->work[p];

				switch (xfrm->op) {
				case XFRM_ADD:
					v += xfrm->args[0];
					break;

				case XFRM_OFST:
					if (v)
						v += xfrm->args[0];
					break;

				case XFRM_ZERO:
					if (v == bottom)
						v = 0.0;
					break;

				case XFRM_MUL: /* negative inverts starting from 1.0 */
					if (xfrm->args[0] >= 0.0)
						v *= xfrm->args[0];
					else
						v = 1.0 + v * xfrm->args[0];
					break;

				case XFRM_GAM:
					if (xfrm->args[0] > 0) {
						if (v < 0.001)
							v = 0;
						else if (v >= 1.000)
							v = 1.0;
						else
							v = pow(v, xfrm->args[0]);
					} else {
						if (v < 0.001)
							v = 0;
						else if (v >= 1.000)
							v = 1.0;
						else
							v = 1.0 - pow(1.0 - v, 1.0 / xfrm->args[0]);
					}
					break;

				case XFRM_GAM2:
					if (v < 0.001)
						v = 0;
					else if (v >= 1.000)
						v = 1.0;
					else if (xfrm->args[0] > 1.0)
						v = fmax(pow(v, xfrm->args[0]), 1.0 - pow(1.0 - v, 1.0 / xfrm->args[0]));
					else if (xfrm->args[0] < 1.0)
						v = fmin(pow(v, xfrm->args[0]), 1.0 - pow(1.0 - v, 1.0 / xfrm->args[0]));
					break;

				case XFRM_HASH:
					if ((x ^ y) & 1)
						v -= soften[p];
					break;

				case XFRM_TWINS:
					/* average two adjacent X pixels, use (0, 2x) for x<128, (2x-255,255) for x>=128 */
					v = soften[p];

					if (v < 0.0)
						v = 0;
					else if (v > 2.0)
						v = 2.0;

					if (v < 1.0) {
						if ((x ^ y) & 1)
							v = v;
						else
							v = 0;
					} else {
						if ((x ^ y) & 1)
							v = 1.0;
						else
							v = v - 1.0;
					}
					if (v < 0.0)
						abort();
					else if (v > 1.0)
						abort();
					break;

				case XFRM_VERTICAL:
					/* average two adjacent X pixels, use (0, 2x) for x<128, (2x-255,255) for x>=128 */
					v = soften[p];

					if (v < 0.0)
						v = 0;
					else if (v > 2.0)
						v = 2.0;

					if (v < 1.0) {
						if (!(x & 1))
							v = v; // strong pixel
						else
							v = 0; // soft pixel
					} else {
						if (!(x & 1))
							v = 1.0; // strong pixel
						else
							v = v - 1.0; // soft pixel
					}
					if (v < 0.0)
						abort();
					else if (v > 1.0)
						abort();
					break;

				case XFRM_QUANTIZE:
					/* make sure to create even steps from 0.0 to 1.0 inclusive */
					if (v > 0)
						v = floor(v * xfrm->args[0]) / (xfrm->args[0] - 1);
					else
						v = 0;
					break;


				case XFRM_QFREQ:
					if (v <= 0.0)
						q = 0;
					else if (v >= 1.0)
						q = 1023;
					else
						q = 1023 * v;

					v = (double)freq[q] / (img->h * img->w);

					/* make sure to create even steps from 0.0 to 1.0 inclusive */
					if (v > 0)
						v = floor(v * xfrm->args[0]) / (xfrm->args[0] - 1);
					else
						v = 0;
					break;

				case XFRM_SOFTEN:
					v -= soften[p] * xfrm->args[0];
					break;

				case XFRM_NORMALIZE:
					v = (v - norm_min) / (norm_max - norm_min);
					break;

				case XFRM_MAP:
					/* map input range to [0..1] */
					if (v < xfrm->args[0])
						v = xfrm->args[0];
					else if (v > xfrm->args[1])
						v = xfrm->args[1];

					v -= xfrm->args[0];
					v /= xfrm->args[1] - xfrm->args[0];

					/* apply gamma2 over the intermediary window */
					if (xfrm->args[4] > 1.0)
						v = fmax(pow(v, xfrm->args[4]), 1.0 - pow(1.0 - v, 1.0 / xfrm->args[4]));
					else if (xfrm->args[4] < 1.0)
						v = fmin(pow(v, xfrm->args[4]), 1.0 - pow(1.0 - v, 1.0 / xfrm->args[4]));

					/* remap to the output range */
					v *= xfrm->args[3] - xfrm->args[2];
					v += xfrm->args[2];
					break;

				default:
					break;
				}

				if (v < 0.0)
					v = 0.0;
				else if (v > 1.0)
					v = 1.0;

				img->work[p] = v;
			}
		}

		if (xfrm->op == XFRM_MAX_LENGTH && (int)xfrm->args[0] > 0) {
			for (y = 0; y < img->h; y++) {
				int line_offset = (int)(((rand() & 255) / 255.0) * xfrm->args[0]);
				float accumulator = 0.0;
				float last = 0.0;
				int length = 0;

				for (x = (y & 1) ? (img->w - 1) : 0;
				     (y & 1) ? (int)x >= 0 : x < img->w;
				     x += (y & 1) ? -1 : 1) {
					uint32_t p = y * img->w + x;
					float v;
					int x2;

					/* only pixels of at least half intensity are limited */
					img->work[p] += last;
					v = img->work[p];
					last = 0;

					if (v < 0.5 || accumulator < 0.85)
						length = 0;
					else
						length++;

					accumulator = ((xfrm->args[0] - 1) * accumulator + v) / xfrm->args[0];

					if (length >= (int)(xfrm->args[0])) {
						/* too many consecutive pixels, remove a previous one */
						x2 = (y & 1) ? x + line_offset : x - line_offset;
						if (x2 < 0)
							x2 = 0;
						else if (x2 >= img->w)
							x2 = img->w - 1;

						last = img->work[y * img->w + x2];
						img->work[y * img->w + x2] = 0;
						length = 0;
					}

				}
			}
		}

		free(soften);
		soften = NULL;
		free(freq);
		freq = NULL;
		xfrm = xfrm->next;
	}
	return 1;
}

/* create a new g-code pass from previous known one */
struct pass *pass_new(struct pass *last)
{
	struct pass *new;

	new = malloc(sizeof(*new));
	if (!new)
		return NULL;

	new->mode    = PASS_MODE_RASTER;
	new->feed    = -1;
	new->spindle = -1;
	new->passes  = 1;
	if (last)
		last->next = new;
	return new;
}

/* emits G-CODE on file <out> (stdout if NULL), for image <img> according to
 * passes descriptions <passes>. Returns non-zero on success, zero on failure.
 * If argc is non-null, a comment line will be emitted recalling the whole
 * command line.
 */
int emit_gcode(const char *out, struct image *img, const struct pass *passes, int argc, char **argv)
{
	FILE *file = stdout;
	const struct pass *pass;
	int pass_num = 0;
	int pass_cnt = 0;
	float pxw = img->mmw / img->w;
	float br = fabs(machine.beam_w / 2.0); // beam radius

	if (out) {
		file = fopen(out, "w");
		if (!file)
			return 0;
	}

	/* must have at least one pass */
	if (!passes) {
		passes = &(const struct pass){
			.mode = PASS_MODE_RASTER,
			.feed = -1,
			.spindle = -1,
			.passes = 1,
		};
	}

	for (pass = passes; pass; pass = pass->next)
		pass_cnt += pass->passes;

	fprintf(file, "; Generated by png2gcode\n");
	fprintf(file, "; Original image size: %ux%u pixels, %.7gx%.7g mm\n", img->w, img->h, f4(img->mmw), f4(img->mmh));
	fprintf(file, "; Work area: %.7g,%.7g -> %.7g,%.7g mm\n", f4(img->orgx), f4(img->orgy), f4(img->orgx+img->mmw), f4(img->orgy+img->mmh));
	fprintf(file, "; Pixel size: %.7gx%.7g mm\n", f4(pxw), f4(img->mmh / img->h));
	if (img->diam > 0)
		fprintf(file, "; Center: 0,0, Max radius: %.7g mm, Diameter: %.7g mm\n", f4(img->diam / 2), f4(img->diam));
	fprintf(file, "; Total of %d passes\n", pass_cnt);
	if (argc > 0) {
		int a, c;
		char *p;

		fprintf(file, "; command line:");
		for (a = 0; a < argc; a++) {
			p = argv[a];
			fputc(' ', file);
			while ((c = *p++)) {
				/* fully escape unprintable characters */
				if (c < ' ' || c == '\'' || c == '\\' || c >= 0x7F)
					fprintf(file, "$'\\x%02x'", c);
				/* quote regular shell delimiters */
				else if ((c >= 0x20 && c <= 0x2a) || c == ';' || c == '?' ||
					 c == '<' || c == '>' || c == '[' || c == ']' ||
					 c == '`' || c == '{' || c == '}' || c == '|' || c == '~')
					fprintf(file, "'%c'", c);
				else
					fputc(c, file);
			}
		}
		fputc('\n', file);
	}

	/* prologue */
	fprintf(file,
		";\n"
		"M5 S0; turn laser off\n"
		"G21; units=mm\n"
		"G90; absolute coordinates\n"
		"G92 X0 Y0 Z0; reset offsets\n"
		"G0 X%.7g Y%.7g; go back to origin\n",
		f4(img->orgx), f4(img->orgy));

	for (pass = passes; pass; pass = pass->next) {
		int loop;
		for (loop = 0; loop < pass->passes; loop++) {
			int base_spindle = pass->spindle;
			int base_feed    = pass->feed;

			pass_num++;

			if (base_feed < 0)
				base_feed = (pass->mode == PASS_MODE_RASTER || pass->mode == PASS_MODE_RASTER_LR || pass->mode == PASS_MODE_CONTOUR) ?
					DEFAULT_FEED_RASTER : DEFAULT_FEED_SHOW;
			if (base_spindle < 0)
				base_spindle = (pass->mode == PASS_MODE_RASTER || pass->mode == PASS_MODE_RASTER_LR || pass->mode == PASS_MODE_CONTOUR) ?
					DEFAULT_SPINDLE_RASTER : DEFAULT_SPINDLE_SHOW;

			fprintf(file, "; pass %d-%d/%d : mode=%s spindle=%d feed=%d\n",
				pass_num, pass_num - loop + pass->passes - 1, pass_cnt,
				pass_mode_names[pass->mode], base_spindle, base_feed);

			if (pass->mode == PASS_MODE_X) {
				fprintf(file, "%s S%d\nG0 X%.7g Y%.7g\nG1 F%d\nG1 X%.7g\nG1 X%.7g\nM5 S0\n",
					machine.laser_on,
					base_spindle, f4(img->orgx), f4(img->orgy),
					base_feed, f4(img->orgx+img->mmw), f4(img->orgx));
			}
			else if (pass->mode == PASS_MODE_Y) {
				fprintf(file, "%s S%d\nG0 X%.7g Y%.7g\nG1 F%d\nG1 Y%.7g\nG1 Y%.7g\nM5 S0\n",
					machine.laser_on,
					base_spindle, f4(img->orgx), f4(img->orgy),
					base_feed, f4(img->orgy+img->mmh), f4(img->orgy));
			}
			else if (pass->mode == PASS_MODE_AXIS) {
				fprintf(file, "%s S%d\nG0 X%.7g Y%.7g\nG1 F%d\nG1 Y%.7g\nG1 Y%.7g\nG1 X%.7g\nG1 X%.7g\nM5 S0\n",
					machine.laser_on,
					base_spindle, f4(img->orgx), f4(img->orgy), base_feed,
					f4(img->orgy+img->mmh), f4(img->orgy),
					f4(img->orgx+img->mmw), f4(img->orgx));
			}
			else if (pass->mode == PASS_MODE_DIAG) {
				fprintf(file, "%s S%d\nG0 X%.7g Y%.7g\nG1 F%d\nG1 X%.7g Y%.7g\nG1 X%.7g Y%.7g\nM5 S0\n",
					machine.laser_on,
					base_spindle, f4(img->orgx), f4(img->orgy), base_feed,
					f4(img->orgx+img->mmw), f4(img->orgy+img->mmh),
					f4(img->orgx), f4(img->orgy));
			}
			else if (pass->mode == PASS_MODE_FRAME) {
				fprintf(file, "%s S%d\nG0 X%.7g Y%.7g\nG1 F%d\n"
					"G1 Y%.7g\n"
					"G1 X%.7g\n"
					"G1 Y%.7g\n"
					"G1 X%.7g\n"
					"M5 S0\n",
					machine.laser_on,
					base_spindle, f4(img->orgx), f4(img->orgy), base_feed,
					f4(img->orgy+img->mmh), f4(img->orgx+img->mmw),
					f4(img->orgy), f4(img->orgx));
			}
			else if (pass->mode == PASS_MODE_CONTOUR) {
				int lx, ly, rx, ry; /* left x,y; right x,y */
				int x, y;
				int cx;

				/* In order to draw the contour, we have 5 steps:
				 *   - spot leftmost, bottom pixel
				 *   - climb on the left
				 *   - spot rightmost, top pixel
				 *   - go down on the right
				 *   - draw to first pixel
				 */

				lx = ly = rx = ry = -1;

				/* spot leftmost dot */
				for (y = 0; ly == -1 && y < img->h; y++) {
					for (x = 0; x < img->w; x++) {
						if (img->work[y * img->w + x]) {
							lx = x;
							ly = y;
							break;
						}
					}
				}

				/* spot rightmost dot */
				for (y = img->h - 1; ry == -1 && y >= 0; y--) {
					for (x = img->w - 1; x >= 0; x--) {
						if (img->work[y * img->w + x]) {
							rx = x;
							ry = y;
							break;
						}
					}
				}

				if (lx < 0 || rx < 0) {
					fprintf(file, "; empty image\n");
					goto contour_end;
				}

				/* OK, start from left */
				fprintf(file, "G0 X%.7g Y%.7g\n%s S%d\nG1 F%d\n",
					f4(img->orgx + imgxr(img, lx)),
					f4(img->orgy + imgyr(img, ly)),
					machine.laser_on,
					base_spindle,
					base_feed);

				cx = lx;
				for (y = ly + 1; y <= ry; y++) {
					for (x = 0; x < img->w; x++) {
						if (img->work[y * img->w + x])
							break;
					}

					if (x == img->w)
						continue;

					if (x != cx)
						fprintf(file, "X%.7g ", f4(img->orgx + imgxr(img, x)));
					fprintf(file, "Y%.7g\n", f4(img->orgy + imgyr(img, y)));
					cx = x;
				}

				/* go top right */
				fprintf(file, "X%.7g Y%.7g\n",
					f4(img->orgx + imgxr(img, rx)),
					f4(img->orgy + imgyr(img, ry)));

				/* go down on the right */
				cx = rx;
				for (y = ry - 1; y >= ly; y--) {
					for (x = img->w - 1; x >= 0; x--) {
						if (img->work[y * img->w + x])
							break;
					}

					if (x < 0)
						continue;

					if (x != cx)
						fprintf(file, "X%.7g ", f4(img->orgx + imgxr(img, x)));
					fprintf(file, "Y%.7g\n", f4(img->orgy + imgyr(img, y)));
					cx = x;
				}

				/* go bottom left */
				fprintf(file, "X%.7g Y%.7g\n",
					f4(img->orgx + imgxr(img, lx)),
					f4(img->orgy + imgyr(img, ly)));

			contour_end:
				fprintf(file, "M5 S0\nG0 X0 Y0\n");
			}
			else if (pass->mode == PASS_MODE_RASTER || pass->mode == PASS_MODE_RASTER_LR) {
				uint32_t curr_spindle;
				unsigned int x, y, x0;
				float xr, yr;   // real positions in millimeters
				float xl; // last non-empty X, in millimeters
				float beam_ofs;
				int ymoved;
				float rl_shift = machine.rl_shift + machine.rl_delay / 1000.0 * pass->feed / 60.0;

				fprintf(file, "%s S%d\nG1 F%d\n",
					machine.laser_on,
					base_spindle, base_feed);

				/* principle: we move even lines from left to right and
				 * odd lines from right to left. Pixels have a width, and
				 * the beam must be lit over a whole pixel. From left to
				 * right, we turn the beam on from the pixel's base position
				 * to the next pixel's base position. From right to left, we
				 * turn the beam on from the next pixel's base position to
				 * the target pixel's base position. The beam is always turned
				 * off when changing line, so each line starts with a G0 which
				 * may be merged with subsequent zeroes. We try to identify
				 * longest lines with same output value and move the beam at
				 * once over them. G0 is used for long zero areas (> 20px).
				 * Note that we don't emit dots but segments. See them as
				 * 1cm-wide segments for an easier representation. Eg:
				 *
				 *    0   1   2   3   4   5   6   7   8   9  <-- px number
				 *  |---|---|---|---|---|---|---|---|---|---|
				 *  0   1   2   3   4   5   6   7   8   9  10 <-- X position
				 *
				 * There is a "from" and a "to" location for each. In L->R
				 * direction, "from" is the previous x (or x0) and "to" is
				 * x. In R->L, the spindle is given by pixel x-1.
				 *
				 * A margin is supported on each side to let the X axis accelerate
				 * and decelerate. This avoids the uneven etching on borders that
				 * particularly affects printed text: overburns under M3, and
				 * underburns under M4. This margin is machine-specific.
				 */
				for (y = 0; y < img->h; y++) {
					/* first pass, left to right */
					if (y >= img->h)
						break;
					yr = y * img->mmh / img->h;
					yr = roundf(yr * 1000.0) / 1000.0;
					ymoved = 1;

					if (br > 0.0 && machine.beam_w < 0.0 && y & 1)
						beam_ofs = br;
					else
						beam_ofs = 0.0;

					curr_spindle = 0;
					x0 = 0;
					xl = -1;
					for (x = 0; x < img->w; x++) {
						uint32_t spindle = img->work[y * img->w + x] * base_spindle;

						xr = x * pxw;
						if (spindle)
							xl = xr;

						if (br > 0.0 && curr_spindle > spindle) {
							/* finish previous pixel to avoid a burn */
							fprintf(file, "X%.7g S%d\n", f4(img->orgx+xr-br+beam_ofs), curr_spindle);
							curr_spindle = spindle;
							if (!spindle)
								x0 = x;
						}

						if (spindle == curr_spindle)
							continue;
						xr = roundf(xr * 1000.0) / 1000.0;
						if (!curr_spindle && (!x0 || ymoved || (x - x0) * pxw >= 2.0 + 2.0 * machine.x_accel)) {
							/* get away at normal speed */
							if (x0 && machine.x_accel)
								fprintf(file, "X%.7g S0\n",
									f4(img->orgx + roundf(x0 * pxw * 1000.0) / 1000.0 + machine.x_accel + br+beam_ofs));

							fprintf(file, "G0 X%.7g", f4(img->orgx+xr+br+beam_ofs-machine.x_accel)); // no lf here, at least one X will follow
							if (ymoved) {
								ymoved = 0;
								fprintf(file, " Y%.7g", f4(img->orgy+yr)); // no lf here, at least one X will follow
							}
							fprintf(file, "\nG1 "); // no lf here, at least one X will follow
							/* finish approach at normal speed */
							if (machine.x_accel)
								fprintf(file, "X%.7g S0\n", f4(img->orgx+xr+br+beam_ofs));
						}
						else
							fprintf(file, "X%.7g S%d\n", f4(img->orgx+xr+br+beam_ofs), curr_spindle);

						curr_spindle = spindle;
						if (!spindle)
							x0 = x;
					}
					/* trace last pixels */
					if (curr_spindle)
						fprintf(file, "X%.7g S%d\n",
							f4(img->orgx - br + beam_ofs + roundf(x * pxw * 1000.0) / 1000.0),
							curr_spindle);

					/* get away at normal speed */
					if (machine.x_accel && xl >= 0)
						fprintf(file, "X%.7g S0\n",
							f4(img->orgx - br + beam_ofs + machine.x_accel + roundf(xl * 1000.0) / 1000.0));

					/* go back to left position if LR mode */
					if (pass->mode == PASS_MODE_RASTER_LR)
						continue;

					/* second pass, right to left */
					y++;
					if (y >= img->h)
						break;
					yr = y * img->mmh / img->h;
					yr = roundf(yr * 1000.0) / 1000.0;
					ymoved = 1;

					if (br > 0.0 && machine.beam_w < 0.0 && y & 1)
						beam_ofs = br;
					else
						beam_ofs = 0.0;

					curr_spindle = 0;
					x0 = 0;
					xl = -1;
					for (x = img->w; x > 0; x--) {
						uint32_t spindle = img->work[y * img->w + x - 1] * base_spindle;

						xr = x * pxw;
						if (spindle)
							xl = xr;

						if (br > 0.0 && curr_spindle > spindle) {
							/* finish previous pixel to avoid a burn */
							fprintf(file, "X%.7g S%d\n", f4(img->orgx + xr + rl_shift + br + beam_ofs), curr_spindle);
							curr_spindle = spindle;
							if (!spindle)
								x0 = x;
						}

						if (spindle == curr_spindle)
							continue;
						xr = roundf(xr * 1000.0) / 1000.0;
						if (!curr_spindle && (!x0 || ymoved || (x0 - x) * pxw >= 2.0 + 2.0 * machine.x_accel)) {
							/* get away at normal speed */
							if (x0 && machine.x_accel)
								fprintf(file, "X%.7g S0\n",
									f4(img->orgx + roundf(x0 * pxw * 1000.0) / 1000.0 - machine.x_accel + br+beam_ofs));

							fprintf(file, "G0 X%.7g", f4(img->orgx + xr + rl_shift - br + beam_ofs + machine.x_accel)); // no lf here, at least one X will follow
							if (ymoved) {
								ymoved = 0;
								fprintf(file, " Y%.7g", f4(img->orgy+yr)); // no lf here, at least one X will follow
							}
							fprintf(file, "\nG1 "); // no lf here, at least one X will follow
							/* finish approach at normal speed */
							if (machine.x_accel)
								fprintf(file, "X%.7g S0\n", f4(img->orgx + xr + rl_shift - br + beam_ofs));
						} else
							fprintf(file, "X%.7g S%d\n", f4(img->orgx + xr - br + beam_ofs + rl_shift), curr_spindle);

						curr_spindle = spindle;
						if (!spindle)
							x0 = x;
					}
					/* trace last pixels */
					if (curr_spindle)
						fprintf(file, "X%.7g S%d\n",
							f4(img->orgx + br + beam_ofs + rl_shift + roundf(x * pxw * 1000.0) / 1000.0),
							curr_spindle);

					/* get away at normal speed */
					if (machine.x_accel && xl >= 0)
						fprintf(file, "X%.7g S0\n",
							f4(img->orgx - br + beam_ofs + rl_shift - machine.x_accel + roundf(xl * 1000.0) / 1000.0));
				}
				// make sure not to draw lines between passes
				fprintf(file, "M5 S0\nG0\n");
			}
		}
	}

	/* epilogue */
	fprintf(file,
		";\n"
		"M5 S0; turn laser off\n"
		"G0 X0 Y0; go back to origin\n"
		);

	return 1;
}

int main(int argc, char **argv)
{
	float pixw = 0, pixh = 0, imgw = 0, imgh = 0, imgd = 0, orgx = 0, orgy = 0;
	int cropx0 = 0, cropy0 = 0, cropx1 = 0, cropy1 = 0, test_sz = 0;
	int ext_l = 0, ext_r = 0, ext_t = 0, ext_b = 0;
	float text_fg_alpha = 1.0, text_fg_value = 1.0;
	float text_bg_alpha = 0, text_bg_value = 0;
	uint32_t arg_auto_crop = 0;
	int arg_raw_preview = 0;
	int text_zoom = 1;
	int txt_w, txt_h;
	enum out_center center_mode = OUT_CNT_NONE;
	enum out_fmt fmt = OUT_FMT_NONE;
	struct pass *curr_pass = NULL;
	struct pass *last_pass = NULL;
	struct pass *passes = NULL;
	struct xfrm *curr = NULL;
	struct xfrm *xfrm = NULL;
	const char  *in  = NULL;
	const char  *out = NULL;
	struct image img;

	while (1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "hi:o:c:f:a:O:g:G:m:q:Q:d:zHtVnrM:S:F:P:w:A:", long_options, &option_index);
		float arg_f = optarg ? atof(optarg) : 0.0;
		int arg_i   = optarg ? atoi(optarg) : 0;

		if (c == -1)
			break;

		switch (c) {
		case 0: /* long option: long_options[option_index] with arg <optarg> */
			break;

		case OPT_AUTO_CROP:
			arg_auto_crop = arg_i;
			break;

		case OPT_CROP_BOTTOM:
			cropy0 = arg_i;
			break;

		case OPT_CROP_LEFT:
			cropx0 = arg_i;
			break;

		case OPT_CROP_RIGHT:
			cropx1 = arg_i;
			break;

		case OPT_CROP_TOP:
			cropy1 = arg_i;
			break;

		case OPT_EXT_BOTTOM:
			ext_b = arg_i;
			break;

		case OPT_EXT_LEFT:
			ext_l = arg_i;
			break;

		case OPT_EXT_RIGHT:
			ext_r = arg_i;
			break;

		case OPT_EXT_TOP:
			ext_t = arg_i;
			break;

		case OPT_TEXT_ZOOM:
			text_zoom = arg_i;
			break;

		case OPT_TEXT_FG_ALPHA:
			text_fg_alpha = arg_f;
			break;

		case OPT_TEXT_FG_VALUE:
			text_fg_value = arg_f;
			break;

		case OPT_TEXT_BG_ALPHA:
			text_bg_alpha = arg_f;
			break;

		case OPT_TEXT_BG_VALUE:
			text_bg_value = arg_f;
			break;

		case OPT_TEXT:
			curr = xfrm_new_text(curr, XFRM_TEXT, optarg, text_zoom, text_fg_alpha, text_fg_value, text_bg_alpha, text_bg_value);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case OPT_SOFTEN:
			curr = xfrm_new(curr, XFRM_SOFTEN, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case OPT_IMGW:
			imgw = arg_f;
			break;

		case OPT_IMGH:
			imgh = arg_f;
			break;

		case OPT_IMGD:
			imgd = arg_f;
			break;

		case OPT_ORGX:
			orgx = arg_f;
			break;

		case OPT_ORGY:
			orgy = arg_f;
			break;

		case OPT_PIXW:
			pixw = arg_f;
			break;

		case OPT_PIXH:
			pixh = arg_f;
			break;

		case OPT_PIXS:
			pixh = pixw = arg_f;
			break;

		case OPT_RL_SHIFT:
			machine.rl_shift = arg_f;
			break;

		case OPT_RL_DELAY:
			machine.rl_delay = arg_f;
			break;

		case OPT_LASER_ON:
			machine.laser_on = optarg;
			break;

		case OPT_TEST:
			test_sz = arg_i;
			break;

		case OPT_MAP: {  /* [from_min:from_max:]to_min:to_max[:gamma2] */
			char *nptr = optarg;
			char *endptr;
			float args[5];
			int narg;

			for (narg = 0; *nptr && narg < 5;) {
				args[narg++] = strtof(nptr, &endptr);
				if (!*endptr)
					break;
				if (*endptr != ':')
					die(1, "unparsable map <%s> (char %d)\n", optarg, endptr - optarg);
				nptr = endptr + 1;
			}

			switch (narg) {
			case 5: /* from_min:from_max:to_min:to_max:gamma2 */
				break;
			case 4: /* from_min:from_max:to_min:to_max */
				args[4] = 1.0; // default gamma2
				break;
			case 3: /* to_min:to_max:gamma2 */
				args[4] = args[2];
				args[3] = args[1];
				args[2] = args[0];
				args[1] = 1.0;
				args[0] = 0.0;
				break;
			case 2: /* to_min:to_max */
				args[4] = 1.0; // default gamma2
				args[3] = args[1];
				args[2] = args[0];
				args[1] = 1.0;
				args[0] = 0.0;
				break;
			default:
				die(1, "missing argument in map <%s>\n", optarg);
			}

			if (args[0] >= args[1])
				die(1, "map <%s>: from_min must be lower than from_max\n", optarg);

			curr = xfrm_new(curr, XFRM_MAP, 5, args);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;
		}

		case OPT_MAX_LENGTH:
			curr = xfrm_new(curr, XFRM_MAX_LENGTH, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'a':
			curr = xfrm_new(curr, XFRM_ADD, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'O':
			curr = xfrm_new(curr, XFRM_OFST, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'g':
			curr = xfrm_new(curr, XFRM_GAM, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'G':
			curr = xfrm_new(curr, XFRM_GAM2, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'H':
			curr = xfrm_new(curr, XFRM_HASH, 0, NULL);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 't':
			curr = xfrm_new(curr, XFRM_TWINS, 0, NULL);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'V':
			curr = xfrm_new(curr, XFRM_VERTICAL, 0, NULL);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'm':
			curr = xfrm_new(curr, XFRM_MUL, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'n':
			curr = xfrm_new(curr, XFRM_NORMALIZE, 0, NULL);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'q':
			if (arg_f <= 1)
				die(1, "quantize levels must be > 1\n");

			curr = xfrm_new(curr, XFRM_QUANTIZE, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'Q':
			if (arg_f <= 1)
				die(1, "quantize levels must be > 1\n");

			curr = xfrm_new(curr, XFRM_QFREQ, 1, &arg_f);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case 'd':
			material.diffusion = arg_f;
			break;

		case 'h':
			usage(0, argv[0]);
			break;

		case 'i' :
			in = optarg;
			break;

		case 'o' :
			out = optarg;
			break;

		case 'c':
			if (strcmp(optarg, "N") == 0)
				center_mode = OUT_CNT_NONE;
			else if (strcmp(optarg, "A") == 0)
				center_mode = OUT_CNT_AXIS;
			else if (strcmp(optarg, "C") == 0)
				center_mode = OUT_CNT_CIRC;
			else
				die(1, "unsupported centering mode %s\n", optarg);
			break;

		case 'f' :
			if (strcmp(optarg, "png") == 0)
				fmt = OUT_FMT_PNG;
			else if (strcmp(optarg, "gcode") == 0)
				fmt = OUT_FMT_GCODE;
			else
				die(1, "unsupported output format %s\n", optarg);
			break;

		case 'M':
			if (!curr_pass) {
				curr_pass = pass_new(last_pass);
				if (!curr_pass)
					die(1, "failed to allocate a new pass\n");
				if (!passes)
					passes = curr_pass;
			}

			for (curr_pass->mode = PASS_MODE_ORIGIN; curr_pass->mode < PASS_MODES; curr_pass->mode++)
				if (strcmp(optarg, pass_mode_names[curr_pass->mode]) == 0)
					break;
			if (curr_pass->mode >= PASS_MODES)
				die(1, "unknown pass mode %s\n", optarg);
			break;

		case 'F':
			if (!curr_pass) {
				curr_pass = pass_new(last_pass);
				if (!curr_pass)
					die(1, "failed to allocate a new pass\n");
				if (!passes)
					passes = curr_pass;
			}
			curr_pass->feed = arg_i;
			break;

		case 'S':
			if (!curr_pass) {
				curr_pass = pass_new(last_pass);
				if (!curr_pass)
					die(1, "failed to allocate a new pass\n");
				if (!passes)
					passes = curr_pass;
			}
			curr_pass->spindle = arg_i;
			break;

		case 'P':
			if (!curr_pass) {
				curr_pass = pass_new(last_pass);
				if (!curr_pass)
					die(1, "failed to allocate a new pass\n");
				if (!passes)
					passes = curr_pass;
			}
			curr_pass->passes = arg_i;
			/* this completes the current pass */
			last_pass = curr_pass;
			curr_pass = NULL;
			break;
		case 'r':
			arg_raw_preview = 1;
			break;
		case 'w':
			machine.beam_w = arg_f;
			break;
		case 'A':
			machine.x_accel = arg_f;
			break;
		case 'z':
			curr = xfrm_new(curr, XFRM_ZERO, 0, NULL);
			if (!curr)
				die(1, "failed to allocate a new transformation\n");
			if (!xfrm)
				xfrm = curr;
			break;

		case ':': /* missing argument */
		case '?': /* unknown option */
			die(1, "");
		}
	}

	if (optind < argc)
		die(1, "unknown argument %s\n", argv[optind]);

	measure_text_size(xfrm, &txt_w, &txt_h);

	if (!in && !test_sz && !((ext_r+ext_l) && (ext_b+ext_t)) && !(txt_w && txt_h)) {
		die(1,
		    "Missing mandatory PNG input file name (-i file), or test pattern (--test),\n"
		    "or text (--text). Alternately you may create an empty area using\n"
		    "--ext-{left,right} and --ext-{bottom,top}.\n"
		    "Use -h for help.\n");
	}

	if (fmt == OUT_FMT_PNG && !out)
		die(1, "missing mandatory PNG output file name (-o file).\nUse -h for help.\n");

	if (fmt == OUT_FMT_NONE) {
		if (!out ||
		    (strlen(out) >= 4 && strcasecmp(out + strlen(out) - 4, ".txt") == 0) ||
		    (strlen(out) >= 6 && strcasecmp(out + strlen(out) - 6, ".gcode") == 0))
			fmt = OUT_FMT_GCODE;
		else if (strlen(out) >= 4 && strcasecmp(out + strlen(out) - 4, ".png") == 0)
			fmt = OUT_FMT_PNG;
		else
			die(1, "missing mandatory output format (-f png ?).\nUse -h for help.\n");
	}

	if ((imgw && pixw) || (imgh && pixh))
		die(1, "not possible to set both image dimensions and pixel dimensions\n");

	if ((imgw || imgh) && imgd)
		die(1, "not possible to set both image dimensions and diameter\n");

	if (test_sz) {
		/* generate a test pattern this size */
		if (!generate_test_pattern(&img, test_sz, test_sz))
			die(2, "failed to generate test pattern.\n");
	} else if (in) {
		if (!read_rgba_file(in, &img))
			die(2, "failed to read file %s\n", in);

		if (arg_auto_crop)
			img.crop_threshold = (uint32_t)255 * arg_auto_crop / 100;

		if (!rgba_to_gray(&img))
			die(3, "failed to convert image to gray\n");
	}

	if (arg_auto_crop) {
		cropx0 = img.minx ? img.minx - 1 : 0;
		cropx1 = (img.maxx < img.w - 1) ? img.w - 2 - img.maxx : 0;
		cropy0 = img.miny ? img.miny - 1 : 0;
		cropy1 = (img.maxy < img.h - 1) ? img.h - 2 - img.maxy : 0;
	}

	if ((cropx0 || cropy0 || cropx1 || cropy1) &&
	    !crop_gray_image(&img, cropx0, cropy0, img.w - 1 - cropx1, img.h - 1 - cropy1))
		die(4, "failed to crop image\n");

	if (ext_l || ext_r || ext_t || ext_b) {
		if (!extend_gray_image(&img, ext_l, ext_r, ext_t, ext_b))
			die(4, "failed to extend image\n");
	}
	else if (!img.w && !img.h && txt_w && txt_h) {
		if (!extend_gray_image(&img, 0, txt_w, txt_h, 0))
			die(4, "failed to extend image to text size\n");
	}

	if (imgw && !imgh)
		imgh = imgw * img.h / img.w;
	else if (!imgw && imgh)
		imgw = imgh * img.w / img.h;

	if (imgw)
		img.mmw = imgw;
	else if (pixw)
		img.mmw = pixw * img.w;
	else
		img.mmw = 0;

	if (imgh)
		img.mmh = imgh;
	else if (pixh)
		img.mmh = pixh * img.h;
	else
		img.mmw = 0;

	img.diam = imgd;

	if (machine.beam_w >= pixw)
		machine.beam_w = 0.9 * pixw;
	else if (machine.beam_w <= -pixw)
		machine.beam_w = -0.9 * pixw;

	if (fmt == OUT_FMT_GCODE && (!img.mmh || !img.mmw) && (center_mode != OUT_CNT_CIRC || !img.diam))
		die(1, "output image width/height/diameter are mandatory in G-CODE output\n");

	if (!gray_to_work(&img))
		die(3, "failed to convert image to work\n");

	img.orgx = orgx;
	img.orgy = orgy;

	if (center_mode == OUT_CNT_AXIS) {
		img.orgx = -img.mmw / 2;
		img.orgy = -img.mmh / 2;
	}
	else if (center_mode == OUT_CNT_CIRC) {
		float ctrx, ctry, rad;

		if (find_center(&img, &ctrx, &ctry, &rad)) {
			img.orgx = -ctrx;
			img.orgy = -ctry;
			img.diam = 2 * rad;
		}
	}

	if (xfrm && !xfrm_apply(&img, xfrm))
		die(6, "failed to apply one transformation to the image\n");

	if (fmt == OUT_FMT_PNG) {
		if (!work_to_gray(&img, arg_raw_preview))
			die(3, "failed to convert image to gray\n");

		if (!write_gray_file(out, &img))
			die(5, "failed to write file %s\n", out);
	}
	else if (fmt == OUT_FMT_GCODE) {
		if (!emit_gcode(out, &img, passes, argc, argv))
			die(6, "failed to emit G-CODE output\n");
	}

	free_image(&img);
	return 0;
}
