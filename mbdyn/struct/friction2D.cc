/* $Header$ */
/* 
 * MBDyn (C) is a multibody analysis code. 
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2023
 *
 * Marco Morandini	<morandini@aero.polimi.it>
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

/* Copyright (C) 2003 Marco Morandini*/

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include <cmath>

#include "mbpar.h"
#include "datamanforward.h"
#include "friction2D.h"
#include "submat.h"

int sign(const doublereal x) {
	if (x >= 0.) {
		return 1;
	} else if (x < 0.) {
		return -1;
	}
	return 0;
};

void
BasicFriction2D::SetValue(DataManager *pDM,
		VectorHandler&X, VectorHandler&XP,
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	NO_OP;
}

ModLugreFriction2D::ModLugreFriction2D(
		const doublereal s0,
		const doublereal s1,
		const doublereal s2,
		const doublereal k,
		const BasicScalarFunction *const ff) : 
sigma0(s0),
sigma1(s1),
sigma2(s2),
kappa(k),
fss(dynamic_cast<const DifferentiableScalarFunction&>(*ff)),
f({0., 0.})
{
	NO_OP;
}

void
ModLugreFriction2D::SetValue(DataManager *pDM,
		VectorHandler&X, 
		VectorHandler&XP,
		SimulationEntity::Hints *ph,
		const unsigned int solution_startdof)
{
	X.PutCoef(solution_startdof+1,f.x[0]/sigma0);
	X.PutCoef(solution_startdof+2,f.x[1]/sigma0);
}

unsigned int ModLugreFriction2D::iGetNumDof(void) const {
	return 1;
};

std::ostream&
ModLugreFriction2D::DescribeDof(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix
		<< "[1]: ModLugreFriction2D state" << std::endl;
}

void
ModLugreFriction2D::DescribeDof(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(1);
	desc[0] = "ModLugreFriction2D state";
}

std::ostream&
ModLugreFriction2D::DescribeEq(std::ostream& out, const char *prefix, bool bInitial) const
{
	return out << prefix
		<< "[1]: ModLugreFriction2D equation" << std::endl;
}

void
ModLugreFriction2D::DescribeEq(std::vector<std::string>& desc, bool bInitial, int i) const
{
	ASSERT(i == -1 || i == 0);
	desc.resize(1);
	desc[0] = "ModLugreFriction2D equation";
}

DofOrder::Order ModLugreFriction2D::GetDofType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction2D::GetDofType");
	return DofOrder::DIFFERENTIAL;
};

DofOrder::Order ModLugreFriction2D::GetEqType(unsigned int i) const {
	ASSERTMSGBREAK(i<iGetNumDof(), "INDEX ERROR in ModLugreFriction2D::GetEqType");
	return DofOrder::DIFFERENTIAL;
};

const doublereal ModLugreFriction2D::fs(const doublereal vm) const {
	return fss(vm); // *sign(v); FIXME
};

const doublereal ModLugreFriction2D::fsd_vm(const doublereal vm) const {
	return fss.ComputeDiff(vm); // *sign(v); FIXME
};

doublereal ModLugreFriction2D::alphatilde(const doublereal zm,
	const doublereal vm) const {

	doublereal zss = fs(vm)/sigma0;
	doublereal zba = kappa*zss;

	if (zm <= zba) {
	} else if ((zba <= zm) && (zm <= zss)) {
		return (0.5*std::sin(M_PI*sigma0*zm/(fs(vm)*(1.-kappa))
			-M_PI*(1.+kappa)/(2*(1.-kappa)))+0.5);
	} else {
		return 1.;
	}
	return 0.;
};

doublereal ModLugreFriction2D::epsilon(const d2D& z,
	const d2D& v) const {

	doublereal eps = (v.x[0]*z.x[0] + v.x[1]*z.x[1] + 1.) / 2.;
	return eps;
};

doublereal ModLugreFriction2D::alpha(const d2D& z,
	const d2D& v) const {
	doublereal zm = std::sqrt(z.x[0]*z.x[0] + z.x[1]*z.x[1]);
	doublereal vm = std::sqrt(v.x[0]*v.x[0] + v.x[1]*v.x[1]);
	return alphatilde(zm, vm) * epsilon(z, v);
}

doublereal ModLugreFriction2D::alphatilded_zm(const doublereal zm,
	const doublereal vm) const {

	doublereal zss = fs(vm)/sigma0;
	doublereal zba = kappa*zss;

	if (zm <= zba) {
		return 0.;
	} else if ((zba <= zm) && (zm <= zss)) {
		return -M_PI/2*std::cos(M_PI*sigma0*zm/fs(vm)/(1.-kappa)
			-M_PI*(1.+kappa)/2./(1.-kappa))
			*sigma0*zm/std::pow(fs(vm),2)/(1.-kappa)
			*fsd_vm(vm);
	} else {
		return 0.;
	}

	return 0.;
};

void ModLugreFriction2D::alphad_v(const d2D& z,
	const d2D& v, ExpandableMatrix& alpha_v) const {

	doublereal zm = std::sqrt(z.x[0]*z.x[0] + z.x[1]*z.x[1]);
	doublereal vm = std::sqrt(v.x[0]*v.x[0] + v.x[1]*v.x[1]);
	doublereal alphat = alphatilde(zm, vm);
	// d2D der;
	// der.x[0] = alphat * z.x[0];
	// der.x[1] = alphat * z.x[1];
	alpha_v.ReDim(1, 1);
	alpha_v.SetBlockDim(1, 2);
	alpha_v.Set(alphat * z.x[0] / 2., 1, 1, 1);
	alpha_v.Set(alphat * z.x[1] / 2., 1, 1, 2);
	return;
};

void ModLugreFriction2D::alphad_z(const d2D& z,
	const d2D& v, ExpandableRowVector& alpha_z) const {

	doublereal zm = std::sqrt(z.x[0]*z.x[0] + z.x[1]*z.x[1]);
	doublereal vm = std::sqrt(v.x[0]*v.x[0] + v.x[1]*v.x[1]);
	doublereal alphat = alphatilde(zm, vm);
	doublereal alphat_dzm = alphatilded_zm(zm, vm);
	doublereal eps = epsilon(z, v);
	// d2D der;
	// der.x[0] = alphat * z.x[0] + eps * alphat_dvm * sign(z.x[0]);
	// der.x[1] = alphat * z.x[1] + eps * alphat_dvm * sign(z.x[1]);
	alpha_z.ReDim(2);
	alpha_z.Set(alphat * v.x[0] / 2. + eps * alphat_dzm * sign(z.x[0]), 1);
	alpha_z.Set(alphat * v.x[1] / 2. + eps * alphat_dzm * sign(z.x[1]), 2);

	return;
};

d2D ModLugreFriction2D::fc(void) const {
	return f;
};

void ModLugreFriction2D::AssRes(
	SubVectorHandler& WorkVec,
	const unsigned int startdof,
	const unsigned int solution_startdof,
	const doublereal F,
	const d2D v,
	const VectorHandler& X,
	const VectorHandler& XP) /*throw(Elem::ChangedEquationStructure)*/ {

	d2D z = {X(solution_startdof+1), X(solution_startdof+2)};
	d2D zp = {XP(solution_startdof+1), XP(solution_startdof+2)};

	doublereal vm = std::sqrt(v.x[0]*v.x[0] + v.x[1]*v.x[1]);
	doublereal fsvm = fs(vm);

	f.x[0] = sigma0*z.x[0] + sigma1*zp.x[0] + sigma2*v.x[0];
	f.x[1] = sigma0*z.x[1] + sigma1*zp.x[1] + sigma2*v.x[1];
	WorkVec.IncCoef(startdof+1,zp.x[0]-v.x[0]+alpha(z,v)*z.x[0]/fsvm*sigma0);
	WorkVec.IncCoef(startdof+2,zp.x[1]-v.x[1]+alpha(z,v)*z.x[1]/fsvm*sigma0);
};

void ModLugreFriction2D::AssJac(
	FullSubMatrixHandler& WorkMat,
	ExpandableMatrix& dfc,
	const unsigned int startdof,
	const unsigned int solution_startdof,
	const doublereal dCoef,
	const doublereal F,
	const d2D v,
	const VectorHandler& X,
	const VectorHandler& XP,
	const ExpandableRowVector& dF,
	const ExpandableMatrix& dv) const {

	d2D z;
	z.x[0] = X(solution_startdof+1);
	z.x[1] = X(solution_startdof+1);
	doublereal vm = std::sqrt(v.x[0]*v.x[0] + v.x[1]*v.x[1]);
	//doublereal zp = XP(solution_startdof+1);
/*
 * 	attrito
 */
 	dfc.ReDim(2, 2);
	dfc.SetBlockDim(1, 2);
	dfc.SetBlockDim(2, 2);
	dfc.SetBlockIdx(1, startdof+1);
	dfc.Set(sigma0*dCoef+sigma1, 1, 1, 1);
	dfc.Set(sigma0*dCoef+sigma1, 2, 1, 2);
	dfc.Set(sigma2, 2, 2, 1);
	dfc.Set(sigma2, 2, 2, 2);
	dfc.Link(2, &dv);
	
/*
 * 	z
 */
 	doublereal alph = alpha(z,v);

	ExpandableRowVector alpha_d_z;
	alphad_z(z, v, alpha_d_z);
	ExpandableMatrix alpha_d_v;
	alphad_v(z, v, alpha_d_v);
	alpha_d_v.Link(1, &dv);

	doublereal fsvm = fs(vm);
	doublereal fsvmdv = fsd_vm(vm);

	// -dot(dz) - alpha * sigma0 / fss * dz
	WorkMat.IncCoef(startdof+1, startdof+1, -1. - alph * sigma0 / fsvm * dCoef);
	WorkMat.IncCoef(startdof+2, startdof+2, -1. - alph * sigma0 / fsvm * dCoef);

	ExpandableMatrix deq_dv;
	deq_dv.ReDim(2, 2);
	deq_dv.SetBlockDim(1, 1);
	deq_dv.SetBlockDim(2, 2);
	// -z * sigma0 / fss * alpha_{/v} * dv
	deq_dv.Set(-z.x[0] * sigma0 / fsvm, 1, 1, 1);
	deq_dv.Set(-z.x[1] * sigma0 / fsvm, 2, 1, 2);
	deq_dv.Link(1, &alpha_d_v);
	// dv + z * alpha * sigma0 / fss^2 * fss_{/v} * dv
	deq_dv.Set(1. + z.x[0] * alph * sigma0 / fsvm / fsvm * fsvmdv * sign(v.x[0]), 1, 2, 1);
	deq_dv.Set(z.x[0] * alph * sigma0 / fsvm / fsvm * fsvmdv * sign(v.x[1]), 1, 2, 2);
	deq_dv.Set(z.x[2] * alph * sigma0 / fsvm / fsvm * fsvmdv * sign(v.x[0]), 2, 2, 1);
	deq_dv.Set(1. + z.x[2] * alph * sigma0 / fsvm / fsvm * fsvmdv * sign(v.x[1]), 2, 2, 2);
	deq_dv.Link(2, &dv);

	deq_dv.Add(WorkMat, startdof+1);

	// -z sigma0/fss * alpha_{/z} * dz
	ExpandableMatrix deq_dz;
	deq_dz.ReDim(2, 1);
	deq_dz.SetBlockDim(1, 2);
	deq_dz.Set(-z.x[0] * sigma0 / fsvm, 1, 1, 1);
	deq_dz.Set(-z.x[1] * sigma0 / fsvm, 2, 1, 2);
	deq_dz.Link(1, &alpha_d_z);
	deq_dz.Add(WorkMat, startdof+1, dCoef);

//	std::cout << alphad_z(z,v) << std::endl;
/*
 * 	callback: dfc[] = df/d{F,v,(z+dCoef,zp)}
 */
};

const OutputHandler::Dimensions
ModLugreFriction2D::GetEquationDimension(integer index) const {
	// DOF == 1
	OutputHandler::Dimensions dimension = OutputHandler::Dimensions::UnknownDimension;

	switch (index)
	{
		case 1:
			dimension = OutputHandler::Dimensions::Velocity;
			break;
	}

	return dimension;
}




//------------------------
d2D SimpleShapeCoefficient2D::Sh_c(void) const {
	return shc;
}
	
d2D SimpleShapeCoefficient2D::Sh_c(
	const d2D f,
	const doublereal F,
	const d2D v) {
	shc = f;
	return shc;
};

void SimpleShapeCoefficient2D::dSh_c(
	ExpandableMatrix& dShc,
	const d2D f,
	const doublereal F,
	const d2D v,
	const ExpandableMatrix& dfc,
	const ExpandableRowVector& dF,
	const ExpandableMatrix& dv) const {
		dShc.ReDim(2, 1);
		dShc.SetBlockDim(1, 2);
		dShc.Set(1., 1, 1, 1);
		dShc.Set(1., 2, 1, 2);
		dShc.Link(1, &dfc);
};



//---------------------------------------

BasicFriction2D *const ParseFriction(MBDynParser& HP,
	DataManager * pDM) 
{
   const char* sKeyWords[] = { 
      "modlugre" "2D",
      NULL
   };
	enum KeyWords { 
	     MODLUGRE2D = 0,
	     LASTKEYWORD
	};
	/* token corrente */
	KeyWords FuncType;
	
	KeyTable K(HP, sKeyWords);
	
	FuncType = KeyWords(HP.IsKeyWord());
	switch (FuncType) {
	case MODLUGRE2D: {
		doublereal sigma0 = HP.GetReal();
		doublereal sigma1 = HP.GetReal();
		doublereal sigma2 = HP.GetReal();
		doublereal kappa = HP.GetReal();
		const BasicScalarFunction*const sf =
			ParseScalarFunction(HP, pDM);
		return new ModLugreFriction2D(sigma0, sigma1, sigma2, kappa, sf);
		break;
	}
	default: {
		silent_cerr("ParseFriction2D(): unrecognized friction type "
				"at line " << HP.GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		break;
	}
	}
	return 0;
};

BasicShapeCoefficient2D *const ParseShapeCoefficient2D(MBDynParser& HP) {
	const char* sKeyWords[] = {
		"simple" "2D",
 		NULL
	};
	enum KeyWords { 
		SIMPLE2D = 0,
		LASTKEYWORD
	};
	/* token corrente */
	KeyWords FuncType;
	
	KeyTable K(HP, sKeyWords);
	
	FuncType = KeyWords(HP.IsKeyWord());
	switch (FuncType) {
	case SIMPLE2D: {
		return new SimpleShapeCoefficient2D();
		break;
	}
	default: {
		silent_cerr("ParseShapeCoefficient2D(): "
			"unrecognized shape coefficient type "
			"at line " << HP.GetLineData() << std::endl);
		throw MBDynParser::ErrGeneric(MBDYN_EXCEPT_ARGS);
		break;
	}
	}
	return 0;
};
