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


// StreamContentCosim - begin

class StreamContentCosim : public StreamContent {
protected:
	const StructNode * pNode;

public:
	StreamContentCosim(const StructNode *pNode,
		StreamContent::Modifier *pMod /* not needed */ );
	virtual ~StreamContentCosim(void);

	void Prepare(void);
	unsigned GetNumChannels(void) const;
};

// #include "geomdata.h"

StreamContentCosim::StreamContentCosim(const StructNode * pNode,
	StreamContent::Modifier *pMod)
: StreamContent(0, pMod), pNode(pNode)
{
	unsigned int size = sizeof(doublereal)*6;

	buf.resize(size);
	memset(&buf[0], 0, size);
	m_pMod->Set(size, &buf[0]);
}

StreamContentCosim::~StreamContentCosim(void)
{
	NO_OP;
}

void
StreamContentCosim::Prepare(void)
{
	doublereal *dbuf = (doublereal *)&buf[0];

	// position
	const Vec3& X = pNode->GetXCurr();

	// TODO
	dbuf[0] = X(1);
	dbuf[1] = X(2);
	dbuf[2] = X(3);

	// orientation
	dbuf += 3*sizeof(doublereal);

	const Mat3x3& R = pNode->GetRCurr();

	dbuf[0] = R(1, 1);
	dbuf[1] = R(1, 2);
	dbuf[2] = R(1, 3);

	m_pMod->Modify();

	ASSERT(curbuf == &buf[buf.size()]);
}

unsigned
StreamContentCosim::GetNumChannels(void) const
{
	return buf.size()/sizeof(doublereal);
}

// StreamContentCosim - end


// Read stream output element content type for cosimulation
struct CosimStreamOutputReader : public StreamOutputContentTypeReader {
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
	if (!SetStreamOutputContentType("cosimulation", rf)) {
		delete rf;
		return -1;
	}

	return 0;
}

