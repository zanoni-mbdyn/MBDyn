/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati	<pierangelo.masarati@polimi.it>
 * Paolo Mantegazza	<paolo.mantegazza@polimi.it>
 *
 * Dipartimento di Ingegneria Aerospaziale - Politecnico di Milano
 * via La Masa, 34 - 20156 Milano, Italy
 * http://www.aero.polimi.it
 *
 * Changing this copyright notice is forbidden.
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation (version 2 of the License).
 * 
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */

#include "mbconfig.h"

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <unistd.h>

extern char *optarg;
extern int optind;

#define BUFFER_SIZE 256 
#define MAX_OUTNAME_LENGTH 73 

/*
@param outname: output file name (up to 72 characters long) 
@param is_modal_displacement: initial modal displacements flag (1: enabled, 0: disabled)
@param is_modal_velocity: initial modal velocities flag (1: enabled, 0: disabled)
@param mass_direction: index of "mass" component to be used (-1: check consistency, 1: x, 2: y, 3: z)  
@return: exit code 
@brief: femgen main routine that convert the mbdyn.tab, mbdyn.mat files into a
		ASCII finite element model file with extension .fem suitable for MBDyn	
		modal joint.  
*/
extern int femgen(const char outname[MAX_OUTNAME_LENGTH], 
					int32_t is_modal_displacement, 
					int32_t is_modal_velocity, 
					int32_t mass_direction);

static void
usage(FILE *outf, int rc)
{
	fprintf(outf,
"usage: femgen [-hd] [-m {cxyz}] [[-o] <outfile>]\n"
"	-d		no initial modal displacements/velocities\n"
"	-h		this message\n"
"	-m <idx>	index of \"mass\" component to be used\n"
"			('x', 'y', 'z'; 'c' to check consistency)\n"
"	-o <outfile>	output file name (up to 72 characters long)\n"
		);
	exit(rc);

}

int main(int argc, char *argv[])
{
	char outname[MAX_OUTNAME_LENGTH] = { ' ' };
	char buf[BUFFER_SIZE];

	int32_t is_modal_displacement = 1;
	int32_t is_modal_velocity = 1; 
	int32_t mass_direction = 0;
	
	for (;;) {
		int opt = getopt(argc, argv, "dhm:o:");
		if (opt == -1) {
			break;
		}

		switch (opt) {
		case 'd':
			is_modal_displacement = 0;
			is_modal_velocity = 0;
			break;

		case '?':
		case 'h':
			usage(stdout, EXIT_SUCCESS);

		case 'm':
			if (optarg[1] != '\0') {
				usage(stderr, EXIT_FAILURE);
			}

			switch (optarg[0]) {
			case 'c':
				mass_direction = -1;
				break;

			case 'x':
				mass_direction = 1;
				break;

			case 'y':
				mass_direction = 2;
				break;

			case 'z':
				mass_direction = 3;
				break;

			default:
				fprintf(stderr, "femgen: unhandled parameter '%c' for option '-m'\n", optarg[0]);
				usage(stderr, EXIT_FAILURE);
			}
			break;

		case 'o': {
			size_t len = strlen(optarg);
			if (len >= sizeof(outname)) {
				fprintf(stderr, "femgen: output file name '%s' too long; trim to 72 bytes or less\n", optarg);
				exit(EXIT_FAILURE);
			}

			strcpy(outname, optarg);
			} break;

		default:
			fprintf(stderr, "femgen: unhandled option '%c'\n", opt);
			usage(stderr, EXIT_FAILURE);
		}
	}

	if (optind < argc) {
		if (outname[0] != ' ') {
			fprintf(stderr, "femgen: output file name already set using '-o' option\n");
			exit(EXIT_FAILURE);

		} else {
			size_t len = strlen(argv[optind]);
			if (len >= sizeof(outname)) {
				fprintf(stderr, "femgen: output file name '%s' too long; trim to 72 bytes or less\n", argv[optind]);
				exit(EXIT_FAILURE);
			}

			strcpy(outname, argv[optind]);
			optind++;
		}
	}

	if (outname[0] == ' ') {
		fprintf(stdout, "Please enter the model name\n");
		fflush(stdout);
		if (fgets(buf, sizeof(buf), stdin) == NULL) {
			fprintf(stderr, "femgen: no model name provided\n");
			exit(EXIT_FAILURE);
		}

		/* trim trailing newline(s) */
		{
			size_t len = strlen(buf);
			while (len > 0 && (buf[len - 1] == '\n' || buf[len - 1] == '\r')) {
				buf[--len] = '\0';
			}
			if (len == 0) {
				fprintf(stderr, "femgen: no model name provided\n");
				exit(EXIT_FAILURE);
			}
			if (len >= sizeof(outname)) {
				fprintf(stderr, "femgen: output file name '%s' too long; trim to 72 bytes or less\n", buf);
				exit(EXIT_FAILURE);
			}
			strcpy(outname, buf);
		}
	}

	if (optind < argc) {
		fprintf(stderr, "femgen: extra args ignored\n");
	}

	return femgen(outname, is_modal_displacement, is_modal_velocity, mass_direction);
}

