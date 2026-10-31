/* Uses the simplified API, thus requires libpng 1.6 or above */
#include <ctype.h>
#include <getopt.h>
#include <math.h>
#include <stdarg.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <png.h>
#include "font5x7.h"

/* default spindle level and feed rate in g-code mode. */
#define DEFAULT_SPINDLE_SHOW     1
#define DEFAULT_SPINDLE_RASTER   255
#define DEFAULT_FEED_SHOW        4800
#define DEFAULT_FEED_RASTER      1200
#define DEFAULT_FEED_CUT         300

/* max number of arguments for a transformation */
#define MAX_XFRM_ARGS 10

/* 5x7 fonts give 6x9 as we keep one pixel around on X and 2 on Y */
#define TEXT_HEIGHT 7
#define TEXT_WIDTH 5
#define TEXT_Y_MARGIN 2
#define TEXT_X_MARGIN 1

/* flags describing each cell of the cut grid */
#define CUT_ARTWORK  0x01 // the cell covers some artwork
#define CUT_KEEP     0x02 // the cell belongs to the piece to cut
#define CUT_OUTSIDE  0x04 // the cell is reachable from the grid's border
#define CUT_SEEN     0x08 // the boundary between this cell and its right one was traced

/* squared distance of cells that are not reachable from any feature cell */
#define CUT_INF 1e20

/* max number of cells of a cut grid using an automatic resolution */
#define CUT_MAX_CELLS (16 * 1024 * 1024)

/* max number of points a fitted line may give up to end where the points stop
 * being exactly aligned, so that it doesn't bend towards the next curve.
 */
#define CUT_LINE_TRIM 4

struct image {
	int32_t w;     // width in pixels
	int32_t h;     // height in pixels
	uint8_t *rgba; // RGBA buffer
	uint8_t *gray; // grayscale buffer
	float *work;   // work area: 0.0 = white, 1.0+ = black
	float orgx;    // image horizontal origin in millimeters
	float orgy;    // image vertical origin in millimeters
	float mmw;     // image width in millimiters
	float mmh;     // image height in millimeters
	float diam;    // max computed diameter (if -cC) or zero
	uint8_t crop_threshold; // crop above this value
	int32_t minx, maxx;
	int32_t miny, maxy;
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
	OPT_UBW,
	OPT_CUT_MARGIN,
	OPT_CUT_SMOOTH,
	OPT_CUT_BRIDGE,
	OPT_CUT_THRESHOLD,
	OPT_CUT_RES,
	OPT_CUT_TOL,
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
	{"raster-ubw",  no_argument,       0, OPT_UBW          },
	{"cut-margin",  required_argument, 0, OPT_CUT_MARGIN   },
	{"cut-smooth",  required_argument, 0, OPT_CUT_SMOOTH   },
	{"cut-bridge",  required_argument, 0, OPT_CUT_BRIDGE   },
	{"cut-threshold", required_argument, 0, OPT_CUT_THRESHOLD },
	{"cut-res",     required_argument, 0, OPT_CUT_RES      },
	{"cut-tol",     required_argument, 0, OPT_CUT_TOL      },
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
	PASS_MODE_CUT,        // cut around the artwork at a distance
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
	[PASS_MODE_CUT]       = "cut",
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

/* G0/G1 mode used in the queue */
enum g_mode {
	G_MODE_INIT = 0,
	G_MODE_G0,
	G_MODE_G1,
};

/* Contains the last programmed move of othe beam (not sent yet). The head is
 * at <from> before the move and will go to <to> in mode <g_mode> with spindle
 * <spindle>.
 */
struct move_queue {
	struct {
		float x, y;
		enum g_mode g_mode;
	} from, to;
	int spindle;
	int queued; // 0=nothing queued
};


/* raster preferences */
int raster_ubw = 0;

/* settings of the cut line around the artwork */
struct cutting {
	float margin;    // min distance between the artwork and the cut line in mm
	float smooth;    // rounding radius in mm, <0 = same as margin
	float bridge;    // width of the bridges joining separate parts in mm, 0 = none, <0 = auto
	float threshold; // pixels above this work value are artwork
	float res;       // cut grid cell size in mm, <=0 = auto
	float tol;       // path simplification tolerance in mm, <0 = auto
} cutting = {
	.margin    = 2.0,
	.smooth    = -1,
	.bridge    = -1,
	.threshold = 0,
	.res       = 0,
	.tol       = -1,
};

/* Grid used to compute the cut line. Cells are squares of <res> millimeters.
 * Cell (x,y) covers image-relative coordinates from x0+x*res to x0+(x+1)*res
 * and from y0+y*res to y0+(y+1)*res. Distances in <dist> are expressed in
 * cells.
 */
struct cut_grid {
	int32_t w;    // width in cells
	int32_t h;    // height in cells
	float res;    // cell size in millimeters
	float x0;     // image-relative X of the grid's left border in millimeters
	float y0;     // image-relative Y of the grid's bottom border in millimeters
	uint8_t *map; // CUT_* flags for each cell
	float *dist;  // squared distances for each cell, see cut_edt()
};

/* closed loop of points, with image-relative coordinates in millimeters. The
 * last point connects back to the first one. Each point starts a segment to
 * the next one, which is straight unless <arc> says otherwise. An arc from a
 * point to itself is a full circle.
 */
struct cut_loop {
	int32_t npts;          // number of points
	int32_t alloc;         // number of allocated points
	float *pts;            // <npts> (x,y) pairs
	uint8_t *arc;          // NULL, or for each point: 0 = straight, 2 = clockwise arc (G2), 3 = counter-clockwise arc (G3)
	float *ctr;            // NULL, or for each point: (x,y) center of the arc starting there
	struct cut_loop *next; // next loop or NULL
};

/* shortest link between two separate parts of a cut grid: cells <p1> and <p2>
 * of parts <part1> and <part2> are <sqdist> squared cells apart.
 */
struct cut_link {
	float sqdist;
	int32_t part1, part2;
	int32_t p1, p2;
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
float get_pix(const struct image *img, int32_t x, int32_t y)
{
	if (y < 0 || y >= img->h)
		return 0;
	if (x < 0 || x >= img->w)
		return 0;
	return img->work[y * img->w + x];
}

/* computes the square of the distance between points 1 and 2 */
float sqdist(int32_t x1, int32_t y1, int32_t x2, int32_t y2, float resx, float resy)
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
	    "  -M --mode    <mode>          pass mode (origin,x,y,axis,diag,frame,raster,raster-lr,cut)\n"
	    "  -F --feed    <value>         feed rate (mm/min, def:4800 frame, 1200 raster, 300 cut)\n"
	    "  -S --spindle <value>         spindle value for intensity 1.0 (def:1 or 255)\n"
	    "  -P --passes  <value>         number of passes with previous parameters (def:1)\n"
	    "Cut line settings (-M cut):\n"
	    "     --cut-margin <mm>         min distance from the artwork to the cut line (def: 2.0)\n"
	    "     --cut-smooth <mm>         rounding radius, fills gaps below twice it (def: margin)\n"
	    "     --cut-bridge <mm>         width of bridges joining parts, 0=none (def: max(3,2*margin))\n"
	    "     --cut-threshold <value>   pixels above this intensity are artwork (def: 0.0)\n"
	    "     --cut-res <mm>            cut computation resolution (def: min(pixel size,0.1))\n"
	    "     --cut-tol <mm>            cut path simplification tolerance (def: res/2)\n"
	    "Material characteristics:\n"
	    "  -d --diffuse <ratio>         adjacent radiation for hash & soften (def:0.18)\n"
	    "Machine settings:\n"
	    "  -w --beam-width <mm>         beam width in millimeters (def: 0, max:90%% pxw, neg=alt)\n"
	    "  -A --x-accel <mm>            X axis acceleration distance in millimeters (def: 0)\n"
	    "     --rl-shift <millimeters>  offset to apply to R->L path in raster mode (def:0)\n"
	    "     --rl-delay <millisec>     time offset to apply to R->L path in raster mode (def:0)\n"
	    "     --laser-on <cmd>          command to turn laser ON (def:M4)\n"
	    "Preferences for raster mode:\n"
	    "     --raster-ubw              apply the beam width correction uniformly (bigger output)\n"
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
int generate_test_pattern(struct image *img, int32_t w, int32_t h)
{
	int32_t x, y, v;

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

/* write the rgba buffer from <img> into file <file>, or to stdout if <file> is
 * NULL, using the same bottom to top direction as write_gray_file(). Returns
 * non-zero on success, otherwise zero.
 */
int write_rgba_file(const char *file, struct image *img)
{
	const int row_stride = -1; // bottom to top
	png_image png_image;
	int ret;

	memset(&png_image, 0, sizeof(png_image));
	png_image.version = PNG_IMAGE_VERSION;
	png_image.width   = img->w;
	png_image.height  = img->h;
	png_image.format  = PNG_FORMAT_RGBA;

	if (file)
		ret = png_image_write_to_file(&png_image, file, 0, img->rgba, row_stride * PNG_IMAGE_ROW_STRIDE(png_image), NULL);
	else
		ret = png_image_write_to_stdio(&png_image, stdout, 0, img->rgba, row_stride * PNG_IMAGE_ROW_STRIDE(png_image), NULL);
	return ret;
}

/* convert image <img> from gray to opaque rgba. The rgba buffer will be
 * allocated and the gray buffer will be freed. Non-zero is returned on
 * success, zero on error.
 */
int gray_to_rgba(struct image *img)
{
	int32_t p;

	if (!img->gray)
		return 0;

	img->rgba = malloc(img->w * img->h * 4);
	if (!img->rgba)
		return 0;

	for (p = 0; p < img->w * img->h; p++) {
		img->rgba[p * 4 + 0] = img->gray[p];
		img->rgba[p * 4 + 1] = img->gray[p];
		img->rgba[p * 4 + 2] = img->gray[p];
		img->rgba[p * 4 + 3] = 255;
	}

	free(img->gray);
	img->gray = NULL;
	return 1;
}

/* convert image <img> from rgba to gray. Transparent is turned to white. The
 * gray buffer will be allocated and the rgba buffer will be freed. Non-zero is
 * returned on success, zero on error.
 */
int rgba_to_gray(struct image *img)
{
	int32_t x, y;
	int32_t minx, maxx, miny, maxy;

	if (!img->rgba)
		return 0;

	img->gray = malloc(img->w * img->h * sizeof(*img->gray));
	if (!img->gray)
		return 0;

	maxx = 0; minx = img->w - 1;
	maxy = 0; miny = img->h - 1;

	for (y = 0; y < img->h; y++) {
		for (x = 0; x < img->w; x++) {
			int32_t v, p;
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
	int32_t x, y;

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
	int32_t x, y;

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
int crop_gray_image(struct image *img, int32_t x0, int32_t y0, int32_t x1, int32_t y1)
{
	int row_pre, row_post;
	uint8_t *src, *dst;
	int32_t x, y;

	if (x0 < 0 || x0 >= img->w || x1 < 0 || x1 >= img->w || x0 > x1)
		return 0;

	if (y0 < 0 || y0 >= img->h || y1 < 0 || y1 >= img->h || y0 > y1)
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
		yinit = lines * zoom * (TEXT_HEIGHT + TEXT_Y_MARGIN) - 1;
	else if (yalign == 0) // center
		yinit = (img->h + lines * zoom * (TEXT_HEIGHT + TEXT_Y_MARGIN) + 1) / 2 - 1;
	else // top
		yinit = img->h - TEXT_Y_MARGIN;

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
		if (yinit < zoom * (TEXT_HEIGHT + TEXT_Y_MARGIN) - 1 || yinit >= img->h)
			goto skip;

		/* we're going to print between line and next-1 */
		width = next - line;
		if (xalign < 0)
			xchar = 0;
		else if (xalign == 0)
			xchar = xinit - width * zoom * (TEXT_WIDTH + TEXT_X_MARGIN) / 2;
		else
			xchar = xinit - width * zoom * (TEXT_WIDTH + TEXT_X_MARGIN);

		while (xchar < 0) {
			xchar += zoom * (TEXT_WIDTH + TEXT_X_MARGIN);
			line++;
		}

		while (line < next) {
			unsigned char c = *line++;

			if (xchar > img->w - zoom * (TEXT_WIDTH + TEXT_X_MARGIN))
				break;

			if (c < 0x20 || c > 0x7f)
				c = 0x20;
			c -= 0x20;

			for (y = 0; y < zoom * (TEXT_HEIGHT + TEXT_Y_MARGIN); y++) {
				if (y >= zoom * TEXT_HEIGHT && !linenum)
					break;
				for (x = 0; x < zoom * (TEXT_WIDTH + TEXT_X_MARGIN); x++) {
					if (x >= zoom * TEXT_WIDTH && line >= next)
						break;
					if (x < zoom * TEXT_WIDTH && y < zoom * TEXT_HEIGHT &&
					    font5x7[(uint8_t)c][y / zoom] & (1 << ((TEXT_WIDTH - 1) - (x / zoom))))
						img->work[(yinit - y) * img->w + xchar + x] =
							img->work[(yinit - y) * img->w + xchar + x] * (1.0 - fg_alpha) + fg_alpha * fg_value;
					else
						img->work[(yinit - y) * img->w + xchar + x] =
							img->work[(yinit - y) * img->w + xchar + x] * (1.0 - bg_alpha) + bg_alpha * bg_value;
				}
			}
			xchar += zoom * (TEXT_WIDTH + TEXT_X_MARGIN);
		}
	skip:
		yinit -= zoom * (TEXT_HEIGHT + TEXT_Y_MARGIN);
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
			h = xfrm->args[2] * xfrm->args[3] * (TEXT_HEIGHT + TEXT_Y_MARGIN); // lines * zoom
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

			w = wmax * xfrm->args[3] * (TEXT_WIDTH + TEXT_X_MARGIN); // cols * zoom
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
	int32_t x, y;
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

/* releases the areas allocated for cut grid <grid> */
void cut_grid_free(struct cut_grid *grid)
{
	free(grid->map);
	free(grid->dist);
	memset(grid, 0, sizeof(*grid));
}

/* Initializes cut grid <grid> to cover image <img> plus <pad> millimeters on
 * each side, using square cells of <res> millimeters, and sets CUT_ARTWORK on
 * the cells covering any pixel whose work value is above <threshold>. Pixel
 * (x,y) is burnt from X=x*pxw to X=(x+1)*pxw along the line Y=y*pxh, so it is
 * considered as covering Y from (y-0.5)*pxh to (y+0.5)*pxh. All cells
 * overlapping a pixel are marked, so the artwork is never under-estimated.
 * Returns non-zero on success, 0 on failure. The grid must be released using
 * cut_grid_free() in both cases.
 */
int cut_grid_init(struct cut_grid *grid, const struct image *img, float res, float pad, float threshold)
{
	float pxw = img->mmw / img->w;
	float pxh = img->mmh / img->h;
	int32_t cx0, cx1, cy0, cy1;
	int32_t cx, cy;
	int32_t x, y;
	double cells;

	memset(grid, 0, sizeof(*grid));
	grid->res = res;
	grid->x0  = -pad;
	grid->y0  = -0.5 * pxh - pad;

	cells = ceil((img->mmw + 2 * pad) / res) * ceil((img->mmh + 2 * pad) / res);
	if (cells > INT32_MAX / 2)
		return 0;

	grid->w = ceil((img->mmw + 2 * pad) / res);
	grid->h = ceil((img->mmh + 2 * pad) / res);

	grid->map  = calloc((size_t)grid->w * grid->h, sizeof(*grid->map));
	grid->dist = malloc((size_t)grid->w * grid->h * sizeof(*grid->dist));
	if (!grid->map || !grid->dist)
		return 0;

	for (y = 0; y < img->h; y++) {
		for (x = 0; x < img->w; x++) {
			if (img->work[y * img->w + x] <= threshold)
				continue;

			/* range of cells overlapping this pixel, tolerating
			 * rounding errors on cell borders.
			 */
			cx0 = floor((x * pxw - grid->x0) / res + 0.001);
			cx1 = ceil(((x + 1) * pxw - grid->x0) / res - 0.001);
			cy0 = floor(((y - 0.5) * pxh - grid->y0) / res + 0.001);
			cy1 = ceil(((y + 0.5) * pxh - grid->y0) / res - 0.001);

			for (cy = cy0; cy < cy1; cy++)
				for (cx = cx0; cx < cx1; cx++)
					grid->map[cy * grid->w + cx] |= CUT_ARTWORK;
		}
	}
	return 1;
}

/* Computes into <d> the 1D squared distance transform of the <n> samples of
 * <f>, taken every <stride> entries, that is, for each q, the min over all p
 * of (q-p)^2 + f[p*stride]. <v> and <z>
 * are work areas of <n> and <n>+1 entries respectively. This uses the lower
 * envelope of parabolas algorithm from P. Felzenszwalb and D. Huttenlocher,
 * "Distance Transforms of Sampled Functions", which runs in linear time: each
 * sample is the base of a parabola, the envelope is built from left to right
 * by eliminating the parabolas it hides (z[] holds the boundaries between the
 * visible ones), then it is evaluated at each position.
 */
void edt_1d(const float *f, int32_t stride, float *d, int32_t *v, double *z, int32_t n)
{
	int32_t k, q;
	double s;

	k = 0;
	v[0] = 0;
	z[0] = -CUT_INF;
	z[1] = CUT_INF;

	for (q = 1; q < n; q++) {
		while (1) {
			s = ((f[q * stride] + (double)q * q) - (f[v[k] * stride] + (double)v[k] * v[k])) / (2.0 * (q - v[k]));
			if (s > z[k] || !k)
				break;
			k--;
		}
		k++;
		v[k] = q;
		z[k] = s;
		z[k + 1] = CUT_INF;
	}

	k = 0;
	for (q = 0; q < n; q++) {
		while (z[k + 1] < q)
			k++;
		d[q] = (double)(q - v[k]) * (q - v[k]) + f[v[k] * stride];
	}
}

/* Turns the <dist> plane of grid <grid>, which must contain 0 for feature
 * cells and CUT_INF for other cells, into the squared euclidean distance in
 * cells from each cell to the nearest feature cell. The 2D transform is made
 * of a 1D transform on each row followed by a 1D transform on each column,
 * each going through a temporary line before replacing the original one.
 * Returns non-zero on success, 0 on failure.
 */
int cut_edt(struct cut_grid *grid)
{
	int32_t n = grid->w > grid->h ? grid->w : grid->h;
	int32_t x, y;
	float *d;
	int32_t *v;
	double *z;
	int ret = 0;

	d = malloc(n * sizeof(*d));
	v = malloc(n * sizeof(*v));
	z = malloc((n + 1) * sizeof(*z));
	if (!d || !v || !z)
		goto out;

	for (y = 0; y < grid->h; y++) {
		edt_1d(grid->dist + y * grid->w, 1, d, v, z, grid->w);
		for (x = 0; x < grid->w; x++)
			grid->dist[y * grid->w + x] = d[x];
	}

	for (x = 0; x < grid->w; x++) {
		edt_1d(grid->dist + x, grid->w, d, v, z, grid->h);
		for (y = 0; y < grid->h; y++)
			grid->dist[y * grid->w + x] = d[y];
	}
	ret = 1;
 out:
	free(d);
	free(v);
	free(z);
	return ret;
}

/* fills <nb> with the indexes of the 4 direct neighbours of cell <p> in grid
 * <grid>, or -1 for those which would be outside of the grid.
 */
void cut_neighbours(const struct cut_grid *grid, int32_t p, int32_t *nb)
{
	int32_t x = p % grid->w;
	int32_t y = p / grid->w;

	nb[0] = x > 0 ? p - 1 : -1;
	nb[1] = x < grid->w - 1 ? p + 1 : -1;
	nb[2] = y > 0 ? p - grid->w : -1;
	nb[3] = y < grid->h - 1 ? p + grid->w : -1;
}

/* returns the squared distance between point (<px>,<py>) and the segment
 * going from (<ax>,<ay>) to (<bx>,<by>).
 */
float seg_sqdist(float px, float py, float ax, float ay, float bx, float by)
{
	float dx = bx - ax;
	float dy = by - ay;
	float len = dx * dx + dy * dy;
	float t = 0;

	if (len > 0) {
		t = ((px - ax) * dx + (py - ay) * dy) / len;
		if (t < 0)
			t = 0;
		else if (t > 1)
			t = 1;
	}
	dx = ax + t * dx - px;
	dy = ay + t * dy - py;
	return dx * dx + dy * dy;
}

/* sets CUT_KEEP on all cells of grid <grid> whose center is within <radius>
 * cells of the segment joining the centers of cells <p1> and <p2>.
 */
void cut_paint_segment(struct cut_grid *grid, int32_t p1, int32_t p2, float radius)
{
	int32_t x1 = p1 % grid->w, y1 = p1 / grid->w;
	int32_t x2 = p2 % grid->w, y2 = p2 / grid->w;
	int32_t xmin, xmax, ymin, ymax;
	int32_t x, y;

	xmin = floor((x1 < x2 ? x1 : x2) - radius);
	xmax = ceil((x1 > x2 ? x1 : x2) + radius);
	ymin = floor((y1 < y2 ? y1 : y2) - radius);
	ymax = ceil((y1 > y2 ? y1 : y2) + radius);

	if (xmin < 0)
		xmin = 0;
	if (xmax > grid->w - 1)
		xmax = grid->w - 1;
	if (ymin < 0)
		ymin = 0;
	if (ymax > grid->h - 1)
		ymax = grid->h - 1;

	for (y = ymin; y <= ymax; y++)
		for (x = xmin; x <= xmax; x++)
			if (seg_sqdist(x, y, x1, y1, x2, y2) <= radius * radius)
				grid->map[y * grid->w + x] |= CUT_KEEP;
}

/* compares two cut links by distance, for qsort() */
int cut_link_cmp(const void *a, const void *b)
{
	const struct cut_link *la = a;
	const struct cut_link *lb = b;

	return la->sqdist < lb->sqdist ? -1 : la->sqdist > lb->sqdist;
}

/* returns the root of part <part> in union-find array <root> */
int32_t cut_part_root(int32_t *root, int32_t part)
{
	while (root[part] != part) {
		root[part] = root[root[part]];
		part = root[part];
	}
	return part;
}

/* Joins the separate parts of the CUT_KEEP area of grid <grid> with bridges
 * of <width> cells so that they form a single piece. Parts are the 4-connected
 * areas, so that two parts only touching by a corner are bridged as well. The
 * bridges follow the minimum spanning tree of the parts, where the cost of a
 * link between two parts is the shortest distance between their boundaries:
 * this joins N parts with N-1 bridges, always the shortest ones possible. The
 * closest boundary cells of each pair of parts are found by brute force, which
 * is fine for the few parts usually met. The tree is built using Kruskal's
 * algorithm: links are considered from the shortest to the longest one, and a
 * link is kept when it joins parts that were not joined yet. Returns non-zero
 * on success, 0 on failure.
 */
int cut_bridge(struct cut_grid *grid, float width)
{
	int32_t cells = grid->w * grid->h;
	int32_t *part = NULL;     // part of each cell or -1
	int32_t *stack = NULL;    // cells to visit while labelling
	int32_t *bstart = NULL;   // index of each part's first entry in <bcell>
	int32_t *bcell = NULL;    // boundary cells sorted by part
	int32_t *root = NULL;     // union-find of the parts
	struct cut_link *links = NULL;
	struct cut_link link;
	int32_t nparts, nlinks, nb[4];
	int32_t p, q, i, j, a, b, n;
	int32_t sp, bnd;
	float dx, dy, d;
	int ret = 0;

	part  = malloc(cells * sizeof(*part));
	stack = malloc(cells * sizeof(*stack));
	if (!part || !stack)
		goto out;

	/* label the parts and count their boundary cells */
	nparts = 0;
	bnd = 0;
	for (p = 0; p < cells; p++)
		part[p] = -1;

	for (p = 0; p < cells; p++) {
		if (!(grid->map[p] & CUT_KEEP) || part[p] >= 0)
			continue;

		sp = 0;
		part[p] = nparts;
		stack[sp++] = p;
		while (sp) {
			q = stack[--sp];
			cut_neighbours(grid, q, nb);
			for (n = 0; n < 4; n++) {
				if (nb[n] < 0 || !(grid->map[nb[n]] & CUT_KEEP)) {
					bnd++;
					continue;
				}
				if (part[nb[n]] >= 0)
					continue;
				part[nb[n]] = nparts;
				stack[sp++] = nb[n];
			}
		}
		nparts++;
	}

	if (nparts < 2) {
		ret = 1;
		goto out;
	}

	/* list the boundary cells of each part, i.e. the ones having at
	 * least one neighbour outside of the part. <bnd> above may count
	 * a cell multiple times, so it's only an upper bound.
	 */
	bstart = calloc(nparts + 1, sizeof(*bstart));
	bcell  = malloc(bnd * sizeof(*bcell));
	root   = malloc(nparts * sizeof(*root));
	links  = malloc((size_t)nparts * (nparts - 1) / 2 * sizeof(*links));
	if (!bstart || !bcell || !root || !links)
		goto out;

	for (p = 0; p < cells; p++) {
		if (part[p] < 0)
			continue;
		cut_neighbours(grid, p, nb);
		for (n = 0; n < 4; n++) {
			if (nb[n] < 0 || part[nb[n]] < 0) {
				bstart[part[p] + 1]++;
				break;
			}
		}
	}

	for (i = 0; i < nparts; i++)
		bstart[i + 1] += bstart[i];

	/* root[] temporarily serves as the fill index of each part */
	for (i = 0; i < nparts; i++)
		root[i] = bstart[i];

	for (p = 0; p < cells; p++) {
		if (part[p] < 0)
			continue;
		cut_neighbours(grid, p, nb);
		for (n = 0; n < 4; n++) {
			if (nb[n] < 0 || part[nb[n]] < 0) {
				bcell[root[part[p]]++] = p;
				break;
			}
		}
	}

	/* find the shortest link between each pair of parts */
	nlinks = 0;
	for (i = 0; i < nparts; i++) {
		for (j = i + 1; j < nparts; j++) {
			link.sqdist = CUT_INF;
			link.part1 = i;
			link.part2 = j;
			link.p1 = link.p2 = -1;
			for (a = bstart[i]; a < bstart[i + 1]; a++) {
				for (b = bstart[j]; b < bstart[j + 1]; b++) {
					dx = bcell[a] % grid->w - bcell[b] % grid->w;
					dx *= dx;
					if (dx >= link.sqdist)
						continue;
					dy = bcell[a] / grid->w - bcell[b] / grid->w;
					d = dx + dy * dy;
					if (d < link.sqdist) {
						link.sqdist = d;
						link.p1 = bcell[a];
						link.p2 = bcell[b];
					}
				}
			}
			links[nlinks++] = link;
		}
	}

	/* build the minimum spanning tree and draw its links */
	qsort(links, nlinks, sizeof(*links), cut_link_cmp);
	for (i = 0; i < nparts; i++)
		root[i] = i;

	for (i = 0; i < nlinks; i++) {
		a = cut_part_root(root, links[i].part1);
		b = cut_part_root(root, links[i].part2);
		if (a == b)
			continue;
		root[a] = b;
		cut_paint_segment(grid, links[i].p1, links[i].p2, width / 2);
	}
	ret = 1;
 out:
	free(part);
	free(stack);
	free(bstart);
	free(bcell);
	free(root);
	free(links);
	return ret;
}

/* Performs a morphological closing of the CUT_KEEP area of grid <grid> with a
 * disk of <radius> cells. This is the area left untouched by a ball of this
 * radius rolling all around it: gaps narrower than 2*radius are filled and
 * inner corners are rounded, while the area never shrinks. It is made of a
 * dilation (cells within <radius> of the area) followed by an erosion (cells
 * farther than <radius> from anything outside of the dilated area). Returns
 * non-zero on success, 0 on failure.
 */
int cut_close(struct cut_grid *grid, float radius)
{
	int32_t cells = grid->w * grid->h;
	float r2 = radius * radius;
	int32_t p;

	for (p = 0; p < cells; p++)
		grid->dist[p] = (grid->map[p] & CUT_KEEP) ? 0 : CUT_INF;

	if (!cut_edt(grid))
		return 0;

	for (p = 0; p < cells; p++)
		grid->dist[p] = (grid->dist[p] <= r2) ? CUT_INF : 0;

	if (!cut_edt(grid))
		return 0;

	for (p = 0; p < cells; p++)
		if (grid->dist[p] > r2)
			grid->map[p] |= CUT_KEEP;
	return 1;
}

/* Sets CUT_OUTSIDE on all cells of grid <grid> without CUT_KEEP that can be
 * reached from the bottom left cell through other such cells using
 * 4-connectivity, then sets CUT_KEEP on all cells not having CUT_OUTSIDE, so
 * that holes in the piece are filled. The bottom left cell must not have
 * CUT_KEEP. Returns non-zero on success, 0 on failure.
 */
int cut_fill_outside(struct cut_grid *grid)
{
	int32_t cells = grid->w * grid->h;
	int32_t *stack;
	int32_t p, sp, n;
	int32_t nb[4];

	stack = malloc(cells * sizeof(*stack));
	if (!stack)
		return 0;

	sp = 0;
	grid->map[0] |= CUT_OUTSIDE;
	stack[sp++] = 0;
	while (sp) {
		p = stack[--sp];
		cut_neighbours(grid, p, nb);
		for (n = 0; n < 4; n++) {
			if (nb[n] < 0 || grid->map[nb[n]] & (CUT_KEEP | CUT_OUTSIDE))
				continue;
			grid->map[nb[n]] |= CUT_OUTSIDE;
			stack[sp++] = nb[n];
		}
	}
	free(stack);

	for (p = 0; p < cells; p++)
		if (!(grid->map[p] & CUT_OUTSIDE))
			grid->map[p] |= CUT_KEEP;
	return 1;
}

/* appends point (<x>,<y>) to loop <loop>. Returns non-zero on success, 0 on
 * failure.
 */
int cut_loop_add(struct cut_loop *loop, float x, float y)
{
	float *pts;

	if (loop->npts >= loop->alloc) {
		loop->alloc = loop->alloc ? loop->alloc * 2 : 256;
		pts = realloc(loop->pts, loop->alloc * 2 * sizeof(*pts));
		if (!pts)
			return 0;
		loop->pts = pts;
	}
	loop->pts[2 * loop->npts] = x;
	loop->pts[2 * loop->npts + 1] = y;
	loop->npts++;
	return 1;
}

/* releases all loops of list <loops> */
void cut_loops_free(struct cut_loop *loops)
{
	struct cut_loop *next;

	while (loops) {
		next = loops->next;
		free(loops->pts);
		free(loops->arc);
		free(loops->ctr);
		free(loops);
		loops = next;
	}
}

/* returns the signed angle in radians swept by the segment starting at point
 * <i> of loop <loop>, positive when counter-clockwise, or 0 if the segment is
 * straight.
 */
double cut_seg_sweep(const struct cut_loop *loop, int32_t i)
{
	int32_t j = (i + 1) % loop->npts;
	double sweep;

	if (!loop->arc || !loop->arc[i])
		return 0;

	if (i == j)
		return loop->arc[i] == 3 ? 2 * M_PI : -2 * M_PI;

	sweep = atan2(loop->pts[2 * j + 1] - loop->ctr[2 * i + 1], loop->pts[2 * j] - loop->ctr[2 * i]) -
		atan2(loop->pts[2 * i + 1] - loop->ctr[2 * i + 1], loop->pts[2 * i] - loop->ctr[2 * i]);

	if (loop->arc[i] == 3 && sweep < 0)
		sweep += 2 * M_PI;
	else if (loop->arc[i] == 2 && sweep > 0)
		sweep -= 2 * M_PI;
	return sweep;
}

/* computes into <x>,<y> the position at fraction <t> (0 to 1) of the segment
 * starting at point <i> of loop <loop>, straight or arc.
 */
void cut_seg_point(const struct cut_loop *loop, int32_t i, double t, float *x, float *y)
{
	int32_t j = (i + 1) % loop->npts;
	double r, a;

	if (!loop->arc || !loop->arc[i]) {
		*x = loop->pts[2 * i] + (loop->pts[2 * j] - loop->pts[2 * i]) * t;
		*y = loop->pts[2 * i + 1] + (loop->pts[2 * j + 1] - loop->pts[2 * i + 1]) * t;
		return;
	}

	r = hypot(loop->pts[2 * i] - loop->ctr[2 * i], loop->pts[2 * i + 1] - loop->ctr[2 * i + 1]);
	a = atan2(loop->pts[2 * i + 1] - loop->ctr[2 * i + 1], loop->pts[2 * i] - loop->ctr[2 * i]);
	a += cut_seg_sweep(loop, i) * t;
	*x = loop->ctr[2 * i] + r * cos(a);
	*y = loop->ctr[2 * i + 1] + r * sin(a);
}

/* computes the extent of the non-empty list of loops <loops> into <minx>,
 * <maxx>, <miny> and <maxy>. Arcs are sampled every few degrees, which is
 * accurate enough for reporting and drawing purposes.
 */
void cut_loops_extent(const struct cut_loop *loops, float *minx, float *maxx, float *miny, float *maxy)
{
	int32_t i, k, steps;
	float x, y;

	*minx = *maxx = loops->pts[0];
	*miny = *maxy = loops->pts[1];
	for (; loops; loops = loops->next) {
		for (i = 0; i < loops->npts; i++) {
			steps = loops->arc && loops->arc[i] ? 128 : 1;
			for (k = 0; k < steps; k++) {
				cut_seg_point(loops, i, (double)k / steps, &x, &y);
				*minx = fmin(*minx, x);
				*maxx = fmax(*maxx, x);
				*miny = fmin(*miny, y);
				*maxy = fmax(*maxy, y);
			}
		}
	}
}

/* rotates the points of loop <loop> so that it starts at its bottom-most
 * point, or the leftmost one if there are several, which eases the alignment
 * of the machine on the cut. Returns non-zero on success, 0 on failure.
 */
int cut_loop_set_start(struct cut_loop *loop)
{
	int32_t start, i;
	float *pts;

	start = 0;
	for (i = 1; i < loop->npts; i++) {
		if (loop->pts[2 * i + 1] < loop->pts[2 * start + 1] ||
		    (loop->pts[2 * i + 1] == loop->pts[2 * start + 1] &&
		     loop->pts[2 * i] < loop->pts[2 * start]))
			start = i;
	}

	if (!start)
		return 1;

	pts = malloc(loop->alloc * 2 * sizeof(*pts));
	if (!pts)
		return 0;

	memcpy(pts, loop->pts + 2 * start, (loop->npts - start) * 2 * sizeof(*pts));
	memcpy(pts + 2 * (loop->npts - start), loop->pts, start * 2 * sizeof(*pts));
	free(loop->pts);
	loop->pts = pts;
	return 1;
}

/* returns non-zero if cell (<x>,<y>) of grid <grid> has CUT_KEEP. Cells outside
 * of the grid don't have it.
 */
int cut_keep(const struct cut_grid *grid, int32_t x, int32_t y)
{
	if (x < 0 || x >= grid->w || y < 0 || y >= grid->h)
		return 0;
	return !!(grid->map[y * grid->w + x] & CUT_KEEP);
}

/* Traces the boundaries of the CUT_KEEP area of grid <grid> using marching
 * squares and prepends them as closed loops to list <loops>. The squares have
 * the centers of 4 adjacent cells as corners, numbered counter-clockwise from
 * the bottom left one, and edge k joins corners k and k+1. Each boundary point
 * is the middle of an edge joining a kept corner and a non-kept one. A loop is
 * followed square by square, always entering a square through an edge going
 * from a kept corner to a non-kept one, and leaving it through the next edge
 * going from a non-kept corner to a kept one. This keeps the piece on the
 * left, so that the outer boundaries are counter-clockwise, and it joins kept
 * cells that only touch by a corner, which means that the area outside of the
 * piece is considered 4-connected. The grid's border cells must not have
 * CUT_KEEP. Returns non-zero on success, 0 on failure.
 */
int cut_trace(struct cut_grid *grid, struct cut_loop **loops)
{
	static const int32_t cdx[4] = { 0, 1, 1, 0 };  // corners of a square
	static const int32_t cdy[4] = { 0, 0, 1, 1 };
	static const int32_t ndx[4] = { 0, 1, 0, -1 }; // square across each edge
	static const int32_t ndy[4] = { -1, 0, 1, 0 };
	struct cut_loop *loop;
	int32_t x, y, sx, sy, k;
	int32_t x0, y0, k0;
	int32_t ax, ay, bx, by;
	int32_t d, m;
	int c[4];

	for (y = 0; y < grid->h; y++) {
		for (x = 0; x < grid->w - 1; x++) {
			if (grid->map[y * grid->w + x] & CUT_SEEN)
				continue;
			if (cut_keep(grid, x, y) == cut_keep(grid, x + 1, y))
				continue;

			loop = calloc(1, sizeof(*loop));
			if (!loop)
				return 0;
			loop->next = *loops;
			*loops = loop;

			/* square entered through the edge between these two cells */
			if (cut_keep(grid, x, y)) {
				sx = x;
				sy = y;
				k = 0;
			} else {
				sx = x;
				sy = y - 1;
				k = 2;
			}

			x0 = sx;
			y0 = sy;
			k0 = k;
			do {
				ax = sx + cdx[k];
				ay = sy + cdy[k];
				bx = sx + cdx[(k + 1) % 4];
				by = sy + cdy[(k + 1) % 4];
				if (ay == by)
					grid->map[ay * grid->w + (ax < bx ? ax : bx)] |= CUT_SEEN;

				if (!cut_loop_add(loop,
						  grid->x0 + ((ax + bx) / 2.0 + 0.5) * grid->res,
						  grid->y0 + ((ay + by) / 2.0 + 0.5) * grid->res))
					return 0;

				for (m = 0; m < 4; m++)
					c[m] = cut_keep(grid, sx + cdx[m], sy + cdy[m]);

				for (d = 1; d < 4; d++) {
					m = (k + d) % 4;
					if (!c[m] && c[(m + 1) % 4])
						break;
				}

				sx += ndx[m];
				sy += ndy[m];
				k = (m + 2) % 4;
			} while (sx != x0 || sy != y0 || k != k0);
		}
	}
	return 1;
}

/* Smooths closed loop <loop> by replacing each point with the weighted average
 * (1,2,1)/4 of itself and its neighbours, <passes> times. The traced points
 * follow cell borders, so they zigzag by up to about one cell around the real
 * shape. This noise would prevent arcs from being recognized within a fine
 * tolerance, while two passes are enough to remove most of it, and only move
 * a smooth curve of radius r by about res^2/r. Returns non-zero on success, 0
 * on failure.
 */
int cut_smooth_loop(struct cut_loop *loop, int passes)
{
	int32_t n = loop->npts;
	int32_t i, h, j;
	float *orig;

	if (n < 3)
		return 1;

	orig = malloc(n * 2 * sizeof(*orig));
	if (!orig)
		return 0;

	while (passes--) {
		memcpy(orig, loop->pts, n * 2 * sizeof(*orig));
		for (i = 0; i < n; i++) {
			h = (i + n - 1) % n;
			j = (i + 1) % n;
			loop->pts[2 * i] = (orig[2 * h] + 2 * orig[2 * i] + orig[2 * j]) / 4;
			loop->pts[2 * i + 1] = (orig[2 * h + 1] + 2 * orig[2 * i + 1] + orig[2 * j + 1]) / 4;
		}
	}
	free(orig);
	return 1;
}

/* returns the address of the (x,y) coordinates of point <i> of loop <loop>,
 * modulo the number of points.
 */
const float *cut_pt(const struct cut_loop *loop, int32_t i)
{
	return loop->pts + 2 * (i % loop->npts);
}

/* returns non-zero if all points of loop <loop> between indexes <a> and <b>
 * (modulo the number of points) are within <tol> of the segment joining them.
 */
int cut_fit_line(const struct cut_loop *loop, int32_t a, int32_t b, float tol)
{
	const float *pa = cut_pt(loop, a);
	const float *pb = cut_pt(loop, b);
	const float *pi;
	int32_t i;

	for (i = a + 1; i < b; i++) {
		pi = cut_pt(loop, i);
		if (seg_sqdist(pi[0], pi[1], pa[0], pa[1], pb[0], pb[1]) > tol * tol)
			return 0;
	}
	return 1;
}

/* Checks if the points of loop <loop> from index <a> to <b> (modulo the number
 * of points) lie on an arc within <tol>. If so, the arc's center is set into
 * <cx>,<cy> and its direction into <dir> (2 = clockwise, 3 = counter-clockwise)
 * and non-zero is returned. The arc passes through both ends, so its center is
 * on the chord's bisector. Its position there is the least squares solution of
 * |p-c|^2 = |a-c|^2 for all points p, which is linear once relative to point
 * <a>. The points must also turn monotonically around the center without
 * completing a turn.
 */
int cut_fit_arc(const struct cut_loop *loop, int32_t a, int32_t b, float tol, float *cx, float *cy, uint8_t *dir)
{
	const float *pa = cut_pt(loop, a);
	const float *pb = cut_pt(loop, b);
	const float *pi;
	double mx, my, nx, ny, len;
	double qx, qy, sab, sbb;
	double ccx, ccy, r, t;
	double ang, prev, d, sweep;
	int32_t i;

	/* half chord and unit normal, relative to point a */
	mx = (pb[0] - pa[0]) / 2;
	my = (pb[1] - pa[1]) / 2;
	len = 2 * hypot(mx, my);
	if (len <= 0 || b - a < 2)
		return 0;
	nx = -2 * my / len;
	ny = 2 * mx / len;

	sab = sbb = 0;
	for (i = a + 1; i < b; i++) {
		pi = cut_pt(loop, i);
		qx = pi[0] - pa[0];
		qy = pi[1] - pa[1];
		sab += (qx * qx + qy * qy - 2 * (mx * qx + my * qy)) * 2 * (nx * qx + ny * qy);
		sbb += 4 * (nx * qx + ny * qy) * (nx * qx + ny * qy);
	}

	if (sbb <= 0)
		return 0;

	t = sab / sbb;
	ccx = mx + t * nx;
	ccy = my + t * ny;
	r = hypot(ccx, ccy);
	if (r <= 2 * tol)
		return 0;

	/* follow the angle of each point around the center */
	prev = atan2(-ccy, -ccx);
	sweep = 0;
	for (i = a + 1; i <= b; i++) {
		pi = cut_pt(loop, i);
		qx = pi[0] - pa[0] - ccx;
		qy = pi[1] - pa[1] - ccy;
		if (i < b && fabs(hypot(qx, qy) - r) > tol)
			return 0;

		ang = atan2(qy, qx);
		d = ang - prev;
		if (d > M_PI)
			d -= 2 * M_PI;
		else if (d <= -M_PI)
			d += 2 * M_PI;
		prev = ang;

		/* tolerate small steps backwards due to the noise */
		if (sweep * d < 0 && fabs(d) > tol / r)
			return 0;
		sweep += d;
	}

	if (fabs(sweep) >= 2 * M_PI - 0.1)
		return 0;

	*cx = pa[0] + ccx;
	*cy = pa[1] + ccy;
	*dir = sweep > 0 ? 3 : 2;
	return 1;
}

/* returns the farthest index b, from <a> + 1 to <end>, such that the points of
 * loop <loop> from <a> to <b> fit a line (<arc> = 0) or an arc (<arc> != 0)
 * within <tol>, or <a> if none does. Aligned points are considered as fitting
 * an arc of infinite radius, so that an arc may start with a straight part.
 * The fit is checked at exponentially growing distances, then the limit is
 * refined by binary search. This assumes that a fitting range remains fitting
 * when shortened, which is close enough to reality and keeps the cost at
 * O(n.log(n)).
 */
int32_t cut_fit_reach(const struct cut_loop *loop, int32_t a, int32_t end, float tol, int arc)
{
	int32_t good, bad, step, b;
	uint8_t dir;
	float cx, cy;

	good = arc ? a : a + 1;
	bad = end + 1;
	for (step = 2; ; step *= 2) {
		b = a + step < end ? a + step : end;
		if (b <= good)
			break;
		if (!cut_fit_line(loop, a, b, tol) &&
		    (!arc || !cut_fit_arc(loop, a, b, tol, &cx, &cy, &dir))) {
			bad = b;
			break;
		}
		good = b;
		if (b == end)
			break;
	}

	while (bad - good > 1) {
		b = (good + bad) / 2;
		if (cut_fit_line(loop, a, b, tol) ||
		    (arc && cut_fit_arc(loop, a, b, tol, &cx, &cy, &dir)))
			good = b;
		else
			bad = b;
	}
	return good;
}

/* Replaces the points of closed loop <loop> with a series of straight segments
 * and arcs, so that no original point is farther than <tol> from them. Starting
 * from the first point, the longest line and the longest arc that fit the next
 * points are searched, and the one reaching farthest is kept, the line being
 * preferred on a tie. An arc is only kept if a line over the same points would
 * deviate by more than twice <tol>, so that nearly straight parts never become
 * arcs of huge radius. A line whose last few points are only accepted thanks
 * to the tolerance ends before them, where they start to leave the line, so
 * that straight parts remain exactly straight instead of being slightly tilted
 * towards the next curve. This repeats from the end of each segment until the
 * loop is closed. The first point is always kept. Returns non-zero on success,
 * 0 on failure.
 */
int cut_fit(struct cut_loop *loop, float tol)
{
	int32_t n = loop->npts;
	int32_t a, bl, ba, bt, m;
	uint8_t *arc;
	float *pts, *ctr;
	float cx, cy;
	uint8_t dir;

	if (n < 4)
		return 1;

	pts = malloc(n * 2 * sizeof(*pts));
	arc = malloc(n * sizeof(*arc));
	ctr = malloc(n * 2 * sizeof(*ctr));
	if (!pts || !arc || !ctr) {
		free(pts);
		free(arc);
		free(ctr);
		return 0;
	}

	for (a = m = 0; a < n; m++) {
		bl = cut_fit_reach(loop, a, n, tol, 0);
		ba = cut_fit_reach(loop, a, n, tol, 1);

		pts[2 * m] = loop->pts[2 * a];
		pts[2 * m + 1] = loop->pts[2 * a + 1];
		arc[m] = 0;
		ctr[2 * m] = ctr[2 * m + 1] = 0;
		if (ba > bl && !cut_fit_line(loop, a, ba, 2 * tol) &&
		    cut_fit_arc(loop, a, ba, tol, &cx, &cy, &dir)) {
			arc[m] = dir;
			ctr[2 * m] = cx;
			ctr[2 * m + 1] = cy;
			a = ba;
		}
		else {
			bt = cut_fit_reach(loop, a, bl, tol / 8, 0);
			a = bl - bt <= CUT_LINE_TRIM ? bt : bl;
		}
	}

	free(loop->pts);
	free(loop->arc);
	free(loop->ctr);
	loop->pts = pts;
	loop->arc = arc;
	loop->ctr = ctr;
	loop->npts = m;
	loop->alloc = n;
	return 1;
}

/* Computes the cut line around the artwork of image <img> according to the
 * <cutting> settings, and returns it in <loops> as a list of closed loops with
 * image-relative coordinates in millimeters. The list is empty when there is
 * no artwork, and contains one loop per separate piece. The steps are:
 *   - mark the cells covering the artwork on a grid larger than the image;
 *   - keep the cells within the margin of the artwork (distance transform);
 *   - join the separate parts with bridges;
 *   - round the area with a ball of the smoothing radius;
 *   - fill the holes so that only outer boundaries remain;
 *   - trace these boundaries, smooth them, make them start at their
 *     bottom-left point and simplify them.
 * Returns non-zero on success, 0 on failure.
 */
int cut_compute(const struct image *img, struct cut_loop **loops)
{
	struct cut_grid grid;
	struct cut_loop *loop;
	float res, pad, smooth, bridge, tol;
	float margin, area;
	int32_t cells, p;
	int ret = 0;

	*loops = NULL;
	smooth = cutting.smooth >= 0 ? cutting.smooth : cutting.margin;
	bridge = cutting.bridge >= 0 ? cutting.bridge : fmax(3.0, 2 * cutting.margin);

	res = cutting.res;
	if (res <= 0) {
		res = fmin(fmin(img->mmw / img->w, img->mmh / img->h), 0.1);
		area = (img->mmw + 2 * (cutting.margin + smooth)) * (img->mmh + 2 * (cutting.margin + smooth));
		if (area / (res * res) > CUT_MAX_CELLS)
			res = sqrt(area / CUT_MAX_CELLS);
	}
	tol = cutting.tol >= 0 ? cutting.tol : res / 2;

	/* the closing must not reach the grid's border */
	pad = cutting.margin + smooth + 2 * res;
	if (!cut_grid_init(&grid, img, res, pad, cutting.threshold))
		goto out;

	cells = grid.w * grid.h;
	for (p = 0; p < cells; p++)
		grid.dist[p] = (grid.map[p] & CUT_ARTWORK) ? 0 : CUT_INF;

	if (!cut_edt(&grid))
		goto out;

	margin = cutting.margin / res;
	for (p = 0; p < cells; p++)
		if (grid.dist[p] <= margin * margin)
			grid.map[p] |= CUT_KEEP;

	if (bridge > 0 && !cut_bridge(&grid, bridge / res))
		goto out;

	if (smooth > 0 && !cut_close(&grid, smooth / res))
		goto out;

	if (!cut_fill_outside(&grid))
		goto out;

	if (!cut_trace(&grid, loops))
		goto out;

	for (loop = *loops; loop; loop = loop->next) {
		if (!cut_smooth_loop(loop, 2))
			goto out;
		if (!cut_loop_set_start(loop))
			goto out;
		if (!cut_fit(loop, tol))
			goto out;
	}
	ret = 1;
 out:
	if (!ret) {
		cut_loops_free(*loops);
		*loops = NULL;
	}
	cut_grid_free(&grid);
	return ret;
}

/* returns non-zero if any pass of list <passes> is a cut pass */
int has_cut_pass(const struct pass *passes)
{
	for (; passes; passes = passes->next)
		if (passes->mode == PASS_MODE_CUT)
			return 1;
	return 0;
}

/* Draws the cut line <loops> in red on image <img>, after extending its gray
 * buffer as needed to fit the line, and turning it to rgba so that the line
 * is distinguished from the artwork. Nothing is changed if there is no loop.
 * Returns non-zero on success, 0 on failure.
 */
int draw_cut_loops(struct image *img, const struct cut_loop *loops)
{
	const struct cut_loop *loop;
	float pxw = img->mmw / img->w;
	float pxh = img->mmh / img->h;
	float minx, maxx, miny, maxy;
	float u0, v0, u1, v1;
	float x0, y0, x1, y1;
	int32_t l, r, t, b;
	int32_t i, k, s, steps, parts;
	int32_t x, y;

	if (!loops)
		return 1;

	cut_loops_extent(loops, &minx, &maxx, &miny, &maxy);

	/* pixel x covers X from x*pxw to (x+1)*pxw, row y covers Y around y*pxh */
	l = minx < 0 ? ceil(-minx / pxw) : 0;
	r = floor(maxx / pxw) >= img->w ? floor(maxx / pxw) - img->w + 1 : 0;
	b = miny < -0.5 * pxh ? ceil(-miny / pxh - 0.5) : 0;
	t = maxy >= (img->h - 0.5) * pxh ? floor(maxy / pxh + 0.5) - img->h + 1 : 0;

	if (!extend_gray_image(img, l, r, t, b))
		return 0;

	if (!gray_to_rgba(img))
		return 0;

	/* arcs are drawn as series of short lines of about one pixel */
	for (loop = loops; loop; loop = loop->next) {
		for (i = 0; i < loop->npts; i++) {
			parts = 1;
			if (loop->arc && loop->arc[i])
				parts = ceil(fabs(cut_seg_sweep(loop, i)) *
					     hypot(loop->pts[2 * i] - loop->ctr[2 * i], loop->pts[2 * i + 1] - loop->ctr[2 * i + 1]) /
					     fmin(pxw, pxh)) + 1;

			for (k = 0; k < parts; k++) {
				cut_seg_point(loop, i, (double)k / parts, &x0, &y0);
				cut_seg_point(loop, i, (double)(k + 1) / parts, &x1, &y1);
				u0 = x0 / pxw + l;
				v0 = y0 / pxh + 0.5 + b;
				u1 = x1 / pxw + l;
				v1 = y1 / pxh + 0.5 + b;
				steps = ceil(fmax(fabs(u1 - u0), fabs(v1 - v0))) + 1;
				for (s = 0; s <= steps; s++) {
					x = floor(u0 + (u1 - u0) * s / steps);
					y = floor(v0 + (v1 - v0) * s / steps);
					if (x < 0 || x >= img->w || y < 0 || y >= img->h)
						continue;
					img->rgba[(y * img->w + x) * 4 + 0] = 255;
					img->rgba[(y * img->w + x) * 4 + 1] = 0;
					img->rgba[(y * img->w + x) * 4 + 2] = 0;
				}
			}
		}
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
	new->next    = NULL;
	if (last)
		last->next = new;
	return new;
}

/* Tries to merge contiguous moves at the same spindle value, or emits a
 * possibly pending one. For G1 moves, the spindle value is also inspected.
 * The expected G mode (G_MODE_G0/G1) must be passed in g_mode. G_MODE_INIT
 * flushes any pending move. Spindle is ignored (assumed 0) for INIT and G1.
 * Returns 0.
 */
int gcode_queue_move(FILE *file, struct move_queue *queue, enum g_mode g_mode, float x, float y, int spindle)
{
	//printf("%d X%.7g Y%.7g S%d\n", ((int)g_mode) - 1, x, y, spindle);

	/* first, if there's something queued at a different G mode or spindle, we
	 * must send it now. Note that we can't have anything queued at G_MODE_INIT.
	 */
	if (g_mode < G_MODE_G1)
		spindle = 0;

	/* If only the spindle changes and no move is made, there's no
	 * engraving, it's just an artefact that must be eliminated.
	 */
	if (queue->queued && g_mode == queue->to.g_mode &&
	    x == queue->to.x && y == queue->to.y)
		return 0;

	/* Let's merge G0 moves regardless of the direction of changes. However
	 * for G1, we only merge contiguous moves, i.e. Y changes when there's
	 * already an X move programmed, or X changes when there's already an Y
	 * move programmed.
	 */
	if (queue->queued &&
	    (spindle != queue->spindle ||
	     g_mode != queue->to.g_mode ||
	     (queue->to.g_mode != G_MODE_G0 &&
	      ((queue->from.x != queue->to.x && y != queue->to.y) ||
	       (queue->from.y != queue->to.y && x != queue->to.x))))) {

		if (queue->to.g_mode != queue->from.g_mode)
			fprintf(file, "%s ", queue->to.g_mode == G_MODE_G0 ? "G0" : "G1");

		if (queue->from.g_mode >= G_MODE_G0 && queue->to.y == queue->from.y)
			fprintf(file, "X%.7g", queue->to.x);
		else if (queue->from.g_mode >= G_MODE_G0 && queue->to.x == queue->from.x)
			fprintf(file, "Y%.7g", queue->to.y);
		else
			fprintf(file, "X%.7g Y%.7g", queue->to.x, queue->to.y);

		if (queue->to.g_mode == G_MODE_G1)
			fprintf(file, " S%d", queue->spindle);

		fprintf(file, "\n");
		queue->queued = 0;
		queue->from = queue->to;
	}

	/* same mode and spindle, let's merge */
	queue->to.g_mode = g_mode;
	queue->spindle = spindle;
	if (g_mode != G_MODE_INIT) {
		queue->to.x = x;
		queue->to.y = y;
		queue->queued++;
	}

	return 0;
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
	struct move_queue move_queue = { 0 };
	struct cut_loop *cut_loops = NULL;
	struct cut_loop *cut;
	char *now_str, *p;
	time_t now;
	int i;

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

	if (has_cut_pass(passes) && !cut_compute(img, &cut_loops))
		return 0;

	now = time(NULL);
	now_str = ctime(&now);
	p = strchr(now_str, '\n');
	if (p)
		*p = 0;

	fprintf(file, "; Generated by png2gcode on %s\n", now_str);
	fprintf(file, "; Original image size: %ux%u pixels, %.7gx%.7g mm\n", img->w, img->h, f4(img->mmw), f4(img->mmh));
	fprintf(file, "; Work area: %.7g,%.7g -> %.7g,%.7g mm\n", f4(img->orgx), f4(img->orgy), f4(img->orgx+img->mmw), f4(img->orgy+img->mmh));
	fprintf(file, "; Pixel size: %.7gx%.7g mm\n", f4(pxw), f4(img->mmh / img->h));
	if (cut_loops) {
		float minx, maxx, miny, maxy;

		cut_loops_extent(cut_loops, &minx, &maxx, &miny, &maxy);
		fprintf(file, "; Cut area: %.7g,%.7g -> %.7g,%.7g mm\n",
			f4(img->orgx + minx), f4(img->orgy + miny), f4(img->orgx + maxx), f4(img->orgy + maxy));
	}
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

			if (base_feed < 0 && pass->mode == PASS_MODE_CUT)
				base_feed = DEFAULT_FEED_CUT;
			else if (base_feed < 0)
				base_feed = (pass->mode == PASS_MODE_RASTER || pass->mode == PASS_MODE_RASTER_LR || pass->mode == PASS_MODE_CONTOUR) ?
					DEFAULT_FEED_RASTER : DEFAULT_FEED_SHOW;
			if (base_spindle < 0)
				base_spindle = (pass->mode == PASS_MODE_RASTER || pass->mode == PASS_MODE_RASTER_LR ||
						pass->mode == PASS_MODE_CONTOUR || pass->mode == PASS_MODE_CUT) ?
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
			else if (pass->mode == PASS_MODE_CUT) {
				/* each loop is cut separately, starting and ending at
				 * its first point. The motion mode (G1/G2/G3) is only
				 * repeated when it changes. Arc centers are given
				 * relative to the arc's start (I/J), computed from the
				 * rounded positions so that the start and end radii
				 * remain identical.
				 */
				if (!cut_loops)
					fprintf(file, "; empty image\n");

				for (cut = cut_loops; cut; cut = cut->next) {
					int g, prev_g = 1;

					fprintf(file, "G0 X%.7g Y%.7g\n%s S%d\nG1 F%d\n",
						f4(img->orgx + cut->pts[0]),
						f4(img->orgy + cut->pts[1]),
						machine.laser_on,
						base_spindle,
						base_feed);

					for (i = 0; i < cut->npts; i++) {
						int j = (i + 1) % cut->npts;

						g = (cut->arc && cut->arc[i]) ? cut->arc[i] : 1;
						if (g != prev_g)
							fprintf(file, "G%d ", g);
						prev_g = g;

						fprintf(file, "X%.7g Y%.7g",
							f4(img->orgx + cut->pts[2 * j]),
							f4(img->orgy + cut->pts[2 * j + 1]));

						if (g != 1)
							fprintf(file, " I%.7g J%.7g",
								f4(f4(img->orgx + cut->ctr[2 * i]) - f4(img->orgx + cut->pts[2 * i])),
								f4(f4(img->orgy + cut->ctr[2 * i + 1]) - f4(img->orgy + cut->pts[2 * i + 1])));
						fprintf(file, "\n");
					}

					fprintf(file, "M5 S0\n");
				}
			}
			else if (pass->mode == PASS_MODE_RASTER || pass->mode == PASS_MODE_RASTER_LR) {
				int x, y, x0;
				float xr, yr;   // real positions in millimeters
				float xl; // last non-empty X, in millimeters
				float beam_ofs;
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
				 * If the machine's beam width is set to a non-null value, the
				 * engraving will only occur starting at the beam's radius (half
				 * of the width) and finishing before the end minus the beam's
				 * radius. In order to limit burns, this is only applied to pixels
				 * that are stronger than their predecessor or successor, unless the
				 * uniform beam width is enabled, in which case every dot is followed
				 * by a zero intensity line to complete the pixel. Example:
				 * Beam width: ((|))         |=center ()=burning area
				 * Engraving: |--********--| *=on, -=off
				 *
				 * A margin is supported on each side to let the X axis accelerate
				 * and decelerate. This avoids the uneven etching on borders that
				 * particularly affects printed text: overburns under M3, and
				 * underburns under M4. This margin is machine-specific.
				 */
				for (y = 0; y < img->h; y++) {
					/* first pass, left to right */
					yr = y * img->mmh / img->h;
					yr = roundf(yr * 1000.0) / 1000.0;

					if (br > 0.0 && machine.beam_w < 0.0 && y & 1)
						beam_ofs = br;
					else
						beam_ofs = 0.0;

					/* In order to travel to areas to be etched, we're using 3 phases over
					 * empty areas (need to know the width):
					 *  - 1) exit:
					 *       if currently moving, move to min(new_pos,last_pos+x_accel)
					 *  - 2) travel:
					 *       if (new_pos - last_pos > 2*x_accel), travel to new_pos-x_accel
					 *  - 3) approach:
					 *       move to new_pos
					 */

					for (x = 0; x < img->w; x++)
						if (img->work[y * img->w + x])
							break;
					if (x >= img->w)
						goto skip_empty_line;

					gcode_queue_move(file, &move_queue, G_MODE_G0,
							 f4(img->orgx + br + beam_ofs - machine.x_accel),
							 f4(img->orgy + yr), 0);

					x0 = -1; // no 0-spindle pixel seen yet
					xl = -machine.x_accel;

					for (x = 0; x < img->w; x++) {
						uint32_t spindle = img->work[y * img->w + x] * base_spindle;

						xr = x * pxw;
						if (!spindle) {
							x0 = x;
							continue;
						}

						xr = roundf(xr * 1000.0) / 1000.0;

						/* phase 1 */
						if (x0 >= 0 && move_queue.to.g_mode > G_MODE_G0) {
							/* first non-empty pixel since xl (might be < 0 if never) */
							gcode_queue_move(file, &move_queue, G_MODE_G1,
									 f4(img->orgx + fmin(xl + machine.x_accel, xr) + br + beam_ofs),
									 f4(img->orgy + yr),
									 0);
						}

						/* phase 2: travel over areas larger than 2mm + 2*acceleration */
						if (move_queue.to.g_mode <= G_MODE_G0 || // was already traveling, fuse this one
						    (xr - xl) > 2.0 + 2.0 * machine.x_accel) {
							gcode_queue_move(file, &move_queue, G_MODE_G0,
									 f4(img->orgx + xr + br + beam_ofs - machine.x_accel),
									 f4(img->orgy + yr), 0);
						}

						/* phase 3: approach at normal speed. This will be fused with the
						 * phase 1 if no G0 was sent at phase 2. We always want to send it
						 * after a hole so as to mark the restart point.
						 */
						if (machine.x_accel || x0 < 0)
							gcode_queue_move(file, &move_queue, G_MODE_G1,
									 f4(img->orgx + xr + br + beam_ofs),
									 f4(img->orgy + yr),
									 0);

						/* we may need to inject some 0-spindle moves for sub-pixels.
						 * it's the same position as for the approach above. In any
						 * case these would be fused.
						 */
						if (br > 0.0 && raster_ubw)
							gcode_queue_move(file, &move_queue, G_MODE_G1,
									 f4(img->orgx + xr + br + beam_ofs),
									 f4(img->orgy + yr),
									 0);

						gcode_queue_move(file, &move_queue, G_MODE_G1,
								 f4(img->orgx + xr + pxw - br + beam_ofs),
								 f4(img->orgy+yr),
								 spindle);

						xl = xr;
						x0 = -1; // reset info as not sending 0-pixels
					}

					/* get away at normal speed */
					if (machine.x_accel && xl >= 0)
						gcode_queue_move(file, &move_queue, G_MODE_G1,
								 f4(img->orgx - br + beam_ofs + machine.x_accel + roundf(xl * 1000.0) / 1000.0),
								 f4(img->orgy + yr),
								 0);
				skip_empty_line:
					/* go back to left position if LR mode */
					if (pass->mode == PASS_MODE_RASTER_LR)
						continue;

					/* second pass, right to left */

					y++;
					if (y >= img->h)
						break;
					yr = y * img->mmh / img->h;
					yr = roundf(yr * 1000.0) / 1000.0;

					if (br > 0.0 && machine.beam_w < 0.0 && y & 1)
						beam_ofs = br;
					else
						beam_ofs = 0.0;

					for (x = img->w; x > 0; x--)
						if (img->work[y * img->w + x - 1])
							break;
					if (!x) /* empty line */
						continue;

					gcode_queue_move(file, &move_queue, G_MODE_G0,
							 f4(img->orgx + img->w * pxw - br + beam_ofs + machine.x_accel + rl_shift),
							 f4(img->orgy + yr), 0);

					x0 = -1; // no 0-spindle pixel seen yet
					xl = ((img->w - 1) * pxw) + machine.x_accel;

					for (x = img->w; x > 0; x--) {
						uint32_t spindle = img->work[y * img->w + x - 1] * base_spindle;

						xr = x * pxw;
						if (!spindle) {
							x0 = x;
							continue;
						}

						xr = roundf(xr * 1000.0) / 1000.0;
						/* phase 1 */
						if (x0 >= 0 && move_queue.to.g_mode > G_MODE_G0) {
							/* first non-empty pixel since xl (might be < 0 if never) */
							gcode_queue_move(file, &move_queue, G_MODE_G1,
									 f4(img->orgx + fmax(xl - machine.x_accel, xr) - br + beam_ofs + rl_shift),
									 f4(img->orgy + yr),
									 0);
						}

						/* phase 2: travel over areas larger than 2mm + 2*acceleration */
						if (move_queue.to.g_mode <= G_MODE_G0 || // was already traveling, fuse this one
						    (xl - xr) > 2.0 + 2.0 * machine.x_accel) {
							gcode_queue_move(file, &move_queue, G_MODE_G0,
									 f4(img->orgx + xr - br + beam_ofs + machine.x_accel + rl_shift),
									 f4(img->orgy + yr), 0);
						}

						/* phase 3: approach at normal speed. This will be fused with the
						 * phase 1 if no G0 was sent at phase 2. We always want to send it
						 * after a hole so as to mark the restart point.
						 */
						if (machine.x_accel || x0 < 0)
							gcode_queue_move(file, &move_queue, G_MODE_G1,
									 f4(img->orgx + xr - br + beam_ofs + rl_shift),
									 f4(img->orgy + yr),
									 0);

						/* we may need to inject some 0-spindle moves for sub-pixels.
						 * it's the same position as for the approach above. In any
						 * case these would be fused.
						 */
						if (br > 0.0 && raster_ubw)
							gcode_queue_move(file, &move_queue, G_MODE_G1,
									 f4(img->orgx + xr - br + beam_ofs + rl_shift),
									 f4(img->orgy + yr),
									 0);

						gcode_queue_move(file, &move_queue, G_MODE_G1,
								 f4(img->orgx + xr - pxw + br + beam_ofs + rl_shift),
								 f4(img->orgy+yr),
								 spindle);

						xl = xr;
						x0 = -1; // reset info as not sending 0-pixels
					}

					/* get away at normal speed */
					if (machine.x_accel && xl >= 0)
						gcode_queue_move(file, &move_queue, G_MODE_G1,
								 f4(img->orgx - br + beam_ofs + rl_shift - machine.x_accel + roundf(xl * 1000.0) / 1000.0),
								 f4(img->orgy + yr),
								 0);
				}

				/* flush pending moves */
				gcode_queue_move(file, &move_queue, G_MODE_INIT, 0, 0, 0);

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

	cut_loops_free(cut_loops);
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
	struct image img = { 0 };

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

		case OPT_UBW:
			raster_ubw = 1;
			break;

		case OPT_CUT_MARGIN:
			cutting.margin = arg_f;
			break;

		case OPT_CUT_SMOOTH:
			cutting.smooth = arg_f;
			break;

		case OPT_CUT_BRIDGE:
			cutting.bridge = arg_f;
			break;

		case OPT_CUT_THRESHOLD:
			cutting.threshold = arg_f;
			break;

		case OPT_CUT_RES:
			cutting.res = arg_f;
			break;

		case OPT_CUT_TOL:
			cutting.tol = arg_f;
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

	if (machine.beam_w > pixw)
		machine.beam_w = 0.9 * pixw;
	else if (machine.beam_w < -pixw)
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
		struct cut_loop *cut_loops = NULL;

		/* the cut line is drawn on the preview if a cut pass is set */
		if (has_cut_pass(passes)) {
			if (!img.mmw || !img.mmh)
				die(1, "image width/height are mandatory to preview the cut line\n");
			if (!cut_compute(&img, &cut_loops))
				die(6, "failed to compute the cut line\n");
		}

		if (!work_to_gray(&img, arg_raw_preview))
			die(3, "failed to convert image to gray\n");

		if (!draw_cut_loops(&img, cut_loops))
			die(3, "failed to draw the cut line\n");
		cut_loops_free(cut_loops);

		if (img.rgba && !write_rgba_file(out, &img))
			die(5, "failed to write file %s\n", out);
		else if (!img.rgba && !write_gray_file(out, &img))
			die(5, "failed to write file %s\n", out);
	}
	else if (fmt == OUT_FMT_GCODE) {
		if (!emit_gcode(out, &img, passes, argc, argv))
			die(6, "failed to emit G-CODE output\n");
	}

	free_image(&img);
	return 0;
}
