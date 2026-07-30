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

#include "hint_impl.h"

Hint *
ParseHint(DataManager *pDM, std::string s)
{
	if (strncasecmp(s.c_str(), "drive{", STRLENOF("drive{")) == 0) {
		s += "drive{";
		
		size_t	len = strlen(s.c_str());

		if (s.c_str()[len - 1] != '}') {
			return 0;
		}

		// char *sDriveStr = new char[len + 1];
		// memcpy(sDriveStr, s, len + 1);
		// sDriveStr[len - 1] = ';';
		std::string sDriveStr(s);
		s += ";";

		return new DriveHint(sDriveStr);

	} else if (strncasecmp(s.c_str(), "drive3{", STRLENOF("drive3{")) == 0) {
		s += "drive3{";
		
		size_t	len = strlen(s.c_str());

		if (s.c_str()[len - 1] != '}') {
			return 0;
		}

		// char *sDriveStr = new char[len + 1];
		// memcpy(sDriveStr, s, len + 1);
		// sDriveStr[len - 1] = ';';
		std::string sDriveStr(s);
		s += ";";

		return new TplDriveHint<Vec3>(sDriveStr);

	} else if (strncasecmp(s.c_str(), "drive6{", STRLENOF("drive6{")) == 0) {
		s += STRLENOF("drive6{");
		
		size_t	len = strlen(s.c_str());

		if (s.c_str()[len - 1] != '}') {
			return 0;
		}

		// char *sDriveStr = new char[len + 1];
		// memcpy(sDriveStr, s, len + 1);
		// sDriveStr[len - 1] = ';';
		std::string sDriveStr(s);
		s += ";";

		return new TplDriveHint<Vec6>(sDriveStr);
	} 

	return 0;
}

/* ParsableHint - start */

ParsableHint::ParsableHint(const std::string s)
: sHint(s)
{
	NO_OP;
}

ParsableHint::~ParsableHint(void)
{
	// if (sHint != 0) {
	// 	SAFEDELETEARR(sHint);
	// }
}

/* DriveHint - start */

DriveHint::DriveHint(const std::string s)
: ParsableHint(s)
{
	NO_OP;
}

DriveHint::~DriveHint(void)
{
	NO_OP;
}

DriveCaller *
DriveHint::pCreateDrive(DataManager *pDM) const
{
	std::istringstream in(sHint);
	InputStream In(in);

	MBDynParser& HP(pDM->GetMBDynParser());
	InputStream& InOrig(HP.GetInputStream());

	HP.PutInputStream(In);
	HP.ExpectArg();
	HP.SetDataManager(pDM);

	DriveCaller *pDC(HP.GetDriveCaller(false));

	HP.PutInputStream(InOrig);

	return pDC;
}

/* DriveHint - end */


/* TplDriveHint - start */

TplDriveHint3 tdh3(0);
TplDriveHint6 tdh6(0);

/* TplDriveHint - end */
