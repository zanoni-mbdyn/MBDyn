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

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

//#include <cstdint>
#include <unistd.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>

#include <fstream>

#include "userelem.h"
#include "socketstream_out_elem.h"
#include "socketstreamdrive.h"
#include "bufmod.h"
#include "drive_.h"

// Read stream output element content type for cosimulation
struct CosimStreamOutputReader : public StreamOutputContentTypeReader{
	virtual StreamContent* Read(DataManager* pDM, MBDynParser& HP);
};

extern "C"
int module_init(const char *module_name, void *pdm, void *php){

#if 0
	DataManager	*pDM = (DataManager *)pdm;
	MBDynParser	*pHP = (MBDynParser *)php;
#endif
	// Stream output element content type for cosimulation
	StreamOutputContentTypeReader *rf = new CosimStreamOutputReader;
	if (!SetStreamOutputContentType("FlightGear", rf)) {
		delete rf;
		return -1;
	}

	return 0;
}

