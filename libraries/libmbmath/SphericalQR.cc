/*
 * MBDyn (C) is a multibody analysis code.
 * http://www.mbdyn.org
 *
 * Copyright (C) 2003-2023
 *
 * This code is a partial merge of HmFe and MBDyn.
 *
 * Pierangelo Masarati  <pierangelo.masarati@polimi.it>
 * Paolo Mantegazza     <paolo.mantegazza@polimi.it>
 * Marco Morandini  <morandini@aero.polimi.it>
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
#include "matvec3.h"

#include <algorithm>
#include <cmath>
/*
 * Compute the square root of the symmetric matrix A, B = sqrtm(A)
 *
 * a[0] = A(1,1)
 * a[1] = A(1,2)
 * a[2] = A[2,2]
 *
 * b[0] = B(1,1)
 * b[1] = B(1,2)
 * b[2] = B[2,2]
 *
 * Formulae taken from https://en.wikipedia.org/wiki/Square_root_of_a_2_by_2_matrix
 *
 */
void sqrtm2x2S(const doublereal *const a, doublereal *const b) {

    doublereal det = a[0] * a[2] - a[1] * a[1];
    doublereal trace = a[0] + a[2];
    doublereal s = std::sqrt(det);
    doublereal t = std::sqrt(trace + 2. * s);
    b[0] = (a[0] + s) / t;
    b[1] = a[1] / t;
    b[2] = (a[2] + s) / t;


    // std::cout << "(1,1): " << b[0]*b[0] + b[1]*b[1] - a[0] << std::endl;
    // std::cout << "(1,2): " << b[0]*b[1] + b[1]*b[2] - a[1] << std::endl;
    // std::cout << "(2,2): " << b[1]*b[1] + b[2]*b[2] - a[2] << std::endl;
    return;
}

/*
 * Compute the iverse of the symmetric matrix A, B = inv(A)
 *
 * a[0] = A(1,1)
 * a[1] = A(1,2)
 * a[2] = A[2,2]
 *
 * b[0] = B(1,1)
 * b[1] = B(1,2)
 * b[2] = B[2,2]
 *
 */
void inv2x2S(const doublereal *const a, doublereal *const b) {
    doublereal det = a[0]*a[2] - a[1]*a[1]; //t2+t4;
    doublereal inv_det = 1.0 / det;

    b[0] = a[2] * inv_det;
    b[1] = -a[1] * inv_det; // t8;
    b[2] = a[0] * inv_det;

    // std::cout << "(1,1): " << b[0]*a[0] + b[1]*a[1] << std::endl;
    // std::cout << "(1,2): " << b[0]*a[1] + b[1]*a[2] << std::endl;
    // std::cout << "(1,2): " << b[1]*a[1] + b[2]*a[2] << std::endl;
    return;
}

/*
 * Sort in place v in ascending order, and return the sorting indices
 *
 * Equivalent to Matlab's [v, idx] = sort(v)
 */
std::vector<integer> sort_vector(Vec3& v) {
    std::vector<integer> idx({1, 2, 3});
    std::stable_sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) {return std::abs(v(i1)) < std::abs(v(i2));});
    Vec3 r;
    r(1) = v(idx[0]);
    r(2) = v(idx[1]);
    r(3) = v(idx[2]);
    v = r;
    return idx;
}


/*
 * Given the vector r, with r.Norm() > 0, compute an orthogonal matrix Q with
 * Q(1,:) aligned with r
 * Q(2,:) and Q(3,:) tangent to the unit sphere
 *
 * if(update) then minimize the rotation around r in ordert to match as mush a possible Qold
 */
void SpericalQR(const Vec3 & r, Mat3x3 &Q, const bool update = false, const Mat3x3& Qold = Eye3) {
    Vec3 q1 = r / r.Norm();
    Vec3 q2 = q1;
    // std::cout << "q2: " << q2 << std::endl;
    std::vector<integer> sort_idx = sort_vector(q2);
    // std::cout << "sort_idx: " << sort_idx[0] << " " << sort_idx[1] << " " << sort_idx[2] << std::endl;

    // std::cout << q1 << "\n\n" << q2 << "\n\n" << sort_idx[0] <<
    //     " " << sort_idx[1] << " " << sort_idx[2] << std::endl;

    Vec3 v2;
    v2(sort_idx[2]) = q2(2);
    v2(sort_idx[1]) = -q2(3);
    v2(sort_idx[0]) = 0.;
    // std::cout << "v2: " << v2 << std::endl;
    v2 = v2 / v2.Norm();
    // std::cout << q1.Dot(v2) << std::endl;
    Vec3 v3 = q1.Cross(v2);
    Q = Mat3x3(q1, v2, v3);
    // std::cout << "\nortocheck: " << Q.MulMT(Q) << std::endl;

    // std::cout << "update: " << update << std::endl;
    if (update) {
        // std::cout << "old: " << Qold << std::endl;
        doublereal c[2][2], cct[3], sqrtc[3], isqrtc[3], u[2][2];
        for (int i = 0; i < 2; i++) {
            for (int j = 0; j < 2; j++) {
                c[i][j] = Qold.GetVec(i+1).Dot(Q.GetVec(j+1));
            }
        }
        cct[0] = c[0][0] * c[0][0] + c[0][1] * c[0][1];
        cct[1] = c[0][0] * c[1][0] + c[0][1] * c[1][1];
        cct[2] = c[1][0] * c[1][0] + c[1][1] * c[1][1];
        sqrtm2x2S(cct, sqrtc);
        inv2x2S(sqrtc, isqrtc);
        u[0][0] = c[0][0] * isqrtc[0] + c[1][0] * isqrtc[1];
        u[1][0] = c[0][0] * isqrtc[1] + c[1][0] * isqrtc[2];
        u[0][1] = c[0][1] * isqrtc[0] + c[1][1] * isqrtc[1];
        u[1][1] = c[0][1] * isqrtc[1] + c[1][1] * isqrtc[2];
        // std::cout << "bef: " << Q << std::endl;
        // std::cout << "bef1: " << Q.GetVec(1) << std::endl;
        // std::cout << "bef2: " << Q.GetVec(2) << std::endl;
        // std::cout << "bef3: " << Q.GetVec(3) << std::endl;
        Vec3 q1 = Q.GetVec(2) * u[0][0] + Q.GetVec(3) * u[1][0];
        Vec3 q2 = Q.GetVec(2) * u[0][1] + Q.GetVec(3) * u[1][1];
        Q.PutVec(2, q1);
        Q.PutVec(3, q2);
        // std::cout << "aft: " << Q << std::endl;
        // std::cout << "aft1: " << Q.GetVec(1) << std::endl;
        // std::cout << "aft2: " << Q.GetVec(2) << std::endl;
        // std::cout << "aft3: " << Q.GetVec(3) << std::endl;
        // std::cout << "\nortocheck: " << Q.MulMT(Q) << std::endl;
    }
}

// int main(void) {
//     Vec3 r(878., 234., 1123.); r = r / r.Norm();
//     Mat3x3 Q;
//     SpericalQR(r, Q);
//     Mat3x3 Qold = Q;
//     r = r + Vec3(0.08, 0.0158, 0.0038);
//     SpericalQR(r, Q, true, Qold);
//
//     return 0;
// }
