/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 1996-2023
 *
 * Pierangelo Masarati  <masarati@aero.polimi.it>
 * Paolo Mantegazza     <mantegazza@aero.polimi.it>
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

/*
 AUTHOR: Reinhard Resch <mbdyn-user@a1.net>
        Copyright (C) 2022(-2023) all rights reserved.

        The copyright of this code is transferred
        to Pierangelo Masarati and Paolo Mantegazza
        for use in the software MBDyn as described
        in the GNU Public License version 2.1
*/

#include <cassert>
#include <iostream>

#include "mbconfig.h"           /* This goes first in every *.c,*.cc file */

#include "myassert.h"
#include "solidshape.h"
#include "demangle.h"

template <typename ElementType, sp_grad::index_type iDim>
bool bCheckShapeFunction()
{
     std::cout << "element type: \"" << mbdyn_demangle(typeid(ElementType)) << "\"\n";

     bool bRes = true;
     using namespace sp_grad;

     SpColVectorA<doublereal, iDim> r;
     SpColVectorA<doublereal, ElementType::iNumNodes> h;

     for (index_type i = 1; i <= ElementType::iNumNodes; ++i) {
          std::cout << "node: " << i << "\n";

          ElementType::NodalPosition(i, r);
          ElementType::ShapeFunction(r, h);

          std::cout << "r = {" << r << "}\n";
          std::cout << "h = {" << h << "}\n";

          for (index_type j = 1; j <= ElementType::iNumNodes; ++j) {
               if (h(j) != (i == j)) {
                    bRes = false;
               }
          }
     }

     return bRes;
}

template <typename ElementType, sp_grad::index_type iDim>
bool bCheckShapeFunctionUPC()
{
     std::cout << "element type: \"" << ElementType::ElementName() << "\"\n";

     bool bRes = true;
     using namespace sp_grad;

     SpColVectorA<doublereal, iDim> r;
     SpColVectorA<doublereal, ElementType::ElemTypeDisplacement::iNumNodes> h;
     SpColVectorA<doublereal, ElementType::ElemTypePressure::iNumNodes> g;

     for (index_type i = 1; i <= ElementType::ElemTypeDisplacement::iNumNodes; ++i) {
          std::cout << "node: " << i << "\n";

          ElementType::ElemTypeDisplacement::NodalPosition(i, r);
          ElementType::ElemTypeDisplacement::ShapeFunction(r, h);
          ElementType::ElemTypePressure::ShapeFunction(r, g);

          std::cout << "r = {" << r << "}\n";
          std::cout << "h = {" << h << "}\n";
          std::cout << "g = {" << g << "}\n";

          for (index_type j = 1; j <= ElementType::ElemTypeDisplacement::iNumNodes; ++j) {
               if (h(j) != (i == j)) {
                    bRes = false;
               }

               if (i <= ElementType::ElemTypePressure::iNumNodes && j <= ElementType::ElemTypePressure::iNumNodes) {
                    if (g(j) != (i == j)) {
                         bRes = false;
                    }
               }
          }
     }

     return bRes;
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ4)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle4, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ8)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle8, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ9)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle9, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionQ8r)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Quadrangle8r, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT6h)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Triangle6h, 2>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH8u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron8u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH8p)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron8p, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron20u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH27u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron27u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20upc)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Hexahedron20upc, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20ur)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Hexahedron20ur, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionH20upcr)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Hexahedron20upcr, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionP6u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Pentahedron6u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionP15u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Pentahedron15u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionP15upc)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Pentahedron15upc, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT4u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Tetrahedron4u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT10u)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunction<Tetrahedron10u, 3>()));
}

MBDYN_TESTSUITE_TEST(solidshapetest, bCheckShapeFunctionT10upc)
{
     MBDYN_TESTSUITE_ASSERT((bCheckShapeFunctionUPC<Tetrahedron10upc, 3>()));
}

int main(int argc, char* argv[])
{
     MBDYN_TESTSUITE_INIT(&argc, argv);

     return MBDYN_RUN_ALL_TESTS();
}
