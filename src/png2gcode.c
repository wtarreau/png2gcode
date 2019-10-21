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

struct image {
	uint32_t w;    // width in pixels
	uint32_t h;    // height in pixels
	uint8_t *rgba; // RGBA buffer
	uint8_t *gray; // grayscale buffer
	float *work;   // work area: 0.0 = white, 1.0+ = black
};

enum out_fmt {
	OUT_FMT_NONE = 0,
	OUT_FMT_PNG,
};

enum opt {
	OPT_CROP_BOTTOM = 256,
	OPT_CROP_LEFT,
	OPT_CROP_RIGHT,
	OPT_CROP_TOP,
	OPT_SOFTEN,
};

const struct option long_options[] = {
	{"help",        no_argument,       0, 'h'              },
	{"in",          required_argument, 0, 'i'              },
	{"out",         required_argument, 0, 'o'              },
	{"fmt",         required_argument, 0, 'f'              },
	{"add",         required_argument, 0, 'a'              },
	{"mul",         required_argument, 0, 'm'              },
	{"gamma",       required_argument, 0, 'g'              },
	{"quantize",    required_argument, 0, 'q'              },
	{"soften",      required_argument, 0, OPT_SOFTEN       },
	{"crop-bottom", required_argument, 0, OPT_CROP_BOTTOM  },
	{"crop-left",   required_argument, 0, OPT_CROP_LEFT    },
	{"crop-right",  required_argument, 0, OPT_CROP_RIGHT   },
	{"crop-top",    required_argument, 0, OPT_CROP_TOP     },
	{0,             0,                 0, 0                }
};

/* describe one transformation to apply to the image */
enum xfrm_op {
	XFRM_NOP = 0,
	XFRM_ADD,
	XFRM_MUL,
	XFRM_GAM,
	XFRM_QUANTIZE,
	XFRM_SOFTEN,
};

struct xfrm {
	enum xfrm_op op;
	float arg;
	struct xfrm *next;
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

void usage(int code, const char *cmd)
{
	die(code,
	    "Usage: %s [options]* [transformations]*\n"
	    "  -h --help                    show this help\n"
	    "  -f --fmt <format>            output format (png)\n"
	    "  -i --in <file>               input PNG file name\n"
	    "  -o --out <file>              output file name\n"
	    "     --crop-bottom <size>      crop this number of pixels from the bottom\n"
	    "     --crop-left   <size>      crop this number of pixels from the left\n"
	    "     --crop-right  <size>      crop this number of pixels from the right\n"
	    "     --crop-top    <size>      crop this number of pixels from the top\n"
	    "Transformations are series of operations applied to the work area:\n"
	    "  -a --add <value>             add <value> [-1..1] to the intensity\n"
	    "  -g --gamma <value>           apply gamma value <value>\n"
	    "  -m --mul <value>             multiply intensity by <value>\n"
	    "  -q --quantize <levels>       quantize to <levels> levels\n"
	    "     --soften <value>          subtract neighbors' average times <value>\n"
	    "", cmd);
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

	if (!img->rgba)
		return 0;

	img->gray = malloc(img->w * img->h * sizeof(*img->gray));
	if (!img->gray)
		return 0;

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
		}
	}

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
 * returned on success, zero on error.
 */
int work_to_gray(struct image *img)
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
 * but the extra size is not released.
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
	return 1;
}

/* append a transformation to list <curr> */
struct xfrm *xfrm_new(struct xfrm *curr, enum xfrm_op op, float arg)
{
	struct xfrm *new;

	new = malloc(sizeof(*new));
	if (!new)
		return NULL;

	new->op    = op;
	new->arg   = arg;
	new->next  = NULL;
	if (curr)
		curr->next = new;
	return new;
}

/* apply transformations starting at <xfrm> to image <img>. Return non-zero on
 * success, zero on failure.
 */
int xfrm_apply(struct image *img, struct xfrm *xfrm)
{
	uint32_t x, y;
	float *soften = NULL;

	while (xfrm) {
		if (xfrm->op == XFRM_SOFTEN) {
			uint32_t px, py;
			soften = malloc(img->h * img->w * sizeof(*soften));
			if (!soften)
				return 0;

			for (y = 0; y < img->h; y++) {
				for (x = 0; x < img->w; x++) {
					float v = 0.0;

					for (py = y - 1; py != y + 2; py++) {
						for (px = x - 1; px != x + 2; px++) {
							/* clip to image border and don't count current point */
							if (py >= img->h ||
							    px >= img->w ||
							    (px == x && py == y))
								continue;

							v += img->work[py * img->w + px];
						}
					}
					soften[y * img->w + x] = v / 8.0;
				}
			}
		}

		for (y = 0; y < img->h; y++) {
			for (x = 0; x < img->w; x++) {
				uint32_t p = y * img->w + x;
				float    v = img->work[p];

				switch (xfrm->op) {
				case XFRM_ADD:
					v += xfrm->arg;
					break;

				case XFRM_MUL:
					v *= xfrm->arg;
					break;

				case XFRM_GAM:
					v = exp(log(v + 1.0) / xfrm->arg) / exp(log(2.0) / xfrm->arg);
					break;

				case XFRM_QUANTIZE:
					if (v > 0)
						v = 1.0 - floor((1.0 - v) * xfrm->arg) / xfrm->arg;
					else
						v = 1.0 /*- (xfrm->arg - 1.0)*/ / xfrm->arg;
					//v = 1.0 - v / xfrm->arg;
					break;

				case XFRM_SOFTEN:
					v -= soften[p] * xfrm->arg;
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

		if (xfrm->op == XFRM_SOFTEN) {
			free(soften);
			soften = NULL;
		}

		xfrm = xfrm->next;
	}
	return 1;
}

int main(int argc, char **argv)
{
	int cropx0 = 0, cropy0 = 0, cropx1 = 0, cropy1 = 0;
	enum out_fmt fmt = OUT_FMT_NONE;
	struct xfrm *curr = NULL;
	struct xfrm *xfrm = NULL;
	const char  *in  = NULL;
	const char  *out = NULL;
	struct image img;

	while (1) {
		int option_index = 0;
		int c = getopt_long(argc, argv, "hi:o:f:a:g:m:q:", long_options, &option_index);

		if (c == -1)
			break;

		switch (c) {
		case 0: /* long option: long_options[option_index] with arg <optarg> */
			break;

		case OPT_CROP_BOTTOM:
			cropy0 = atoi(optarg);
			break;

		case OPT_CROP_LEFT:
			cropx0 = atoi(optarg);
			break;

		case OPT_CROP_RIGHT:
			cropx1 = atoi(optarg);
			break;

		case OPT_CROP_TOP:
			cropy1 = atoi(optarg);
			break;

		case OPT_SOFTEN:
			curr = xfrm_new(curr, XFRM_SOFTEN, atof(optarg));
			if (!curr)
				die(1, "failed to allocate a new transformation\n", optarg);
			if (!xfrm)
				xfrm = curr;
			break;

		case 'a':
			curr = xfrm_new(curr, XFRM_ADD, atof(optarg));
			if (!curr)
				die(1, "failed to allocate a new transformation\n", optarg);
			if (!xfrm)
				xfrm = curr;
			break;

		case 'g':
			curr = xfrm_new(curr, XFRM_GAM, atof(optarg));
			if (!curr)
				die(1, "failed to allocate a new transformation\n", optarg);
			if (!xfrm)
				xfrm = curr;
			break;

		case 'm':
			curr = xfrm_new(curr, XFRM_MUL, atof(optarg));
			if (!curr)
				die(1, "failed to allocate a new transformation\n", optarg);
			if (!xfrm)
				xfrm = curr;
			break;

		case 'q':
			if (atof(optarg) < 1)
				die(1, "quantize levels must be >= 1\n");

			curr = xfrm_new(curr, XFRM_QUANTIZE, atof(optarg));
			if (!curr)
				die(1, "failed to allocate a new transformation\n", optarg);
			if (!xfrm)
				xfrm = curr;
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

		case 'f' :
			if (strcmp(optarg, "png") == 0)
				fmt = OUT_FMT_PNG;
			else
				die(1, "unsupported output format %s\n", optarg);
			break;
		}
	}

	if (optind < argc)
		die(1, "unknown argument %s\n", argv[optind]);

	if (!in)
		die(1, "missing mandatory PNG input file name (-i file)\n");

	if (fmt == OUT_FMT_PNG && !out)
		die(1, "missing mandatory PNG output file name (-o file)\n");

	if (out && fmt == OUT_FMT_NONE)
		die(1, "missing mandatory output format (-f png ?)\n");

	if (!read_rgba_file(in, &img))
		die(2, "failed to read file %s\n", in);

	if (!rgba_to_gray(&img))
		die(3, "failed to convert image to gray\n");

	if ((cropx0 || cropy0 || cropx1 || cropy1) &&
	    !crop_gray_image(&img, cropx0, cropy0, img.w - 1 - cropx1, img.h - 1 - cropy1))
		die(4, "failed to crop image\n");

	if (!gray_to_work(&img))
		die(3, "failed to convert image to work\n");

	if (xfrm && !xfrm_apply(&img, xfrm))
		die(6, "failed to apply one transformation to the image\n");

	if (fmt == OUT_FMT_PNG) {
		if (!work_to_gray(&img))
			die(3, "failed to convert image to gray\n");

		if (!write_gray_file(out, &img))
			die(5, "failed to write file %s\n", out);
	}

	free_image(&img);
	return 0;
}
