#include <inttypes.h>
#include <stdio.h>
#include "font5x7.h"

int main(int argc, char **argv)
{
	const char *p;
	int x, y;
	char c;

	if (argc > 1) {
		p = argv[1];
		while ((c = *p++)) {
			if (c < 0x20 || c > 0x7f)
				c = 0x20;
			c -= 0x20;

			for (y = 0; y < 7; y++) {
				for (x = 0; x < 5; x++)
					if (font5x7[(uint8_t)c][y] & (1 << (4-x)))
						putchar('#');
					else
						putchar(' ');
				putchar('\n');
			}
			putchar('\n');
		}
	} else {
		for (c = 0; c < 96; c += 8) {
			for (y = 0; y < 7; y++) {
				do {
					for (x = 0; x < 5; x++)
					  if (font5x7[(uint8_t)c][y] & (1 << (4-x)))
							putchar('#');
						else
							putchar(' ');
					putchar(' ');
				} while ((c = ((c & 7) == 7) ? (c & ~7) : (c + 1)) & 7);
				putchar('\n');
			}
			putchar('\n');
		}
	}
	return 0;
}
