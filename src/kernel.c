/* Copyright (C) 2013 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA
 * 02111-1307  USA
 */

#include <R.h>
#include <Rinternals.h>

#include <math.h>

// naive implementation of the kernel function
////////////////////////////////////////////////////////////////////////////////

SEXP exponential_kernel_1d(SEXP x, SEXP y, SEXP l, SEXP var)
{
        R_len_t i, j;
        R_len_t nx = length(x);
        R_len_t ny = length(y);
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rvar = REAL(var);
        double *rans, norm;
        SEXP ans, dim;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 1) {
                error("x has invalid dimension");
        }

        dim = getAttrib(y, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 1) {
                error("y has invalid dimension");
        }

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(var) != 1) {
                error("var is not a scalar");
        }

        /* compute kernel */
        PROTECT(ans = allocMatrix(REALSXP, nx, ny));
        rans = REAL(ans);
        for(i = 0; i < nx; i++) {
                for(j = 0; j < ny; j++) {
                        norm = (rx[i] - ry[j])*(rx[i] - ry[j]);
                        rans[i + nx*j] =
                                (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*norm);
                }
        }
        UNPROTECT(1);

        return(ans);
}

SEXP exponential_kernel_2d(SEXP x, SEXP y, SEXP l, SEXP var)
{
        R_len_t i, j;
        R_len_t nx;
        R_len_t ny;
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rvar = REAL(var);
        double *rans, norm;
        SEXP ans, dim;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("x has invalid dimension");
        }
        nx = INTEGER(dim)[0];

        dim = getAttrib(y, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("y has invalid dimension");
        }
        ny = INTEGER(dim)[0];

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(var) != 1) {
                error("var is not a scalar");
        }

        /* compute kernel */
        PROTECT(ans = allocMatrix(REALSXP, nx, ny));
        rans = REAL(ans);
        for(i = 0; i < nx; i++) {
                for(j = 0; j < ny; j++) {
                        norm = (rx[i + nx*0] - ry[j + ny*0])*(rx[i + nx*0] - ry[j + ny*0])
                             + (rx[i + nx*1] - ry[j + ny*1])*(rx[i + nx*1] - ry[j + ny*1]);
                        rans[i + nx*j] =
                                (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*norm);
                }
        }
        UNPROTECT(1);

        return(ans);
}

// kernel function on spherical coordinate system (e.g. http://en.wikipedia.org/wiki/Spherical_coordinate_system) 
// on a unit circle (r=1);
// theta: azimuthal angle, phi: polar angle; all angles in radians
// direction of fixation: phi = pi/2, theta = pi, i.e. the zenith is pointing upwards
SEXP exponential_kernel_spherical(SEXP phi, SEXP theta, SEXP l, SEXP var)
{
        R_len_t i, j;
        R_len_t nphi;
        R_len_t ntheta;
        double *rphi	= REAL(phi);
        double *rtheta	= REAL(theta);
        double *rl	= REAL(l);
        double *rvar	= REAL(var);
        double *rans, d;
        SEXP ans, dim;

        /* check input */
        dim = getAttrib(phi, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("phi has invalid dimension");
        }
        nphi = INTEGER(dim)[0];

        dim = getAttrib(theta, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("theta has invalid dimension");
        }
        ntheta = INTEGER(dim)[0];

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(var) != 1) {
                error("var is not a scalar");
        }

        /* compute kernel */
        PROTECT(ans = allocMatrix(REALSXP, nphi, ntheta));
        rans = REAL(ans);
        for(i = 0; i < nphi; i++) {
                for(j = 0; j < ntheta; j++) {
			// d: segment on the unit circle
			d = acos( cos(rtheta[i] - rtheta[i+ntheta])*sin(rphi[i])*sin(rphi[i+nphi]) 
				+ cos(rphi[i])*cos(rphi[i+nphi]) );
                        rans[i + nx*j] =
                                (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*d*d);
                }
        }
        UNPROTECT(1);

        return(ans);
}

// kernel functions as sparse symmetric band matrices
////////////////////////////////////////////////////////////////////////////////

SEXP exponential_kernel_1d_sparse(SEXP x, SEXP y, SEXP l, SEXP var, SEXP n, SEXP m)
{
        R_len_t p, q;
        R_len_t nx   = length(x);
        R_len_t ny   = length(y);
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rvar = REAL(var);
        double *rn   = REAL(n);
        double *rm   = REAL(m);
        double *rans, norm;
        double *rtmp;
        SEXP ans, dim, tmp;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 1) {
                error("x has invalid dimension");
        }

        dim = getAttrib(y, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 1) {
                error("y has invalid dimension");
        }

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(var) != 1) {
                error("var is not a scalar");
        }
        if (length(n) != 1) {
                error("n is not a scalar");
        }
        if (length(m) != 1) {
                error("m is not a scalar");
        }

        /* define a list that contains the diagonals */
        PROTECT(ans = allocVector(VECSXP, *rn + *rm + 1));
        /* loop over sub-diagonals */
        for(p = *rn; p > 0; p--) {
                /* define a vector that contains a diagonal */
                PROTECT(tmp = allocVector(REALSXP, fmin(nx-p, ny)));
                rtmp = REAL(tmp);
                /* loop through the diagonal */
                for(q = 0; q < fmin(nx-p, ny); q++) {
                        /* compute norm */
                        norm = (rx[p+q] - ry[q])*(rx[p+q] - ry[q]);
                        /* compute exponential kernel function */
                        rtmp[q] = (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*norm);
                }
                /* insert the diagonal into the list */
                SET_VECTOR_ELT(ans, *rn - p, tmp);
                UNPROTECT(1);
        }
        /* loop over main diagonal and super-diagonals */
        for(q = 0; q <= *rm; q++) {
                /* define a vector that contains a diagonal */
                PROTECT(tmp = allocVector(REALSXP, fmin(nx, ny-q)));
                rtmp = REAL(tmp);
                /* loop through the diagonal */
                for(p = 0; p < fmin(nx, ny-q); p++) {
                        /* compute norm */
                        norm = (rx[p] - ry[p+q])*(rx[p] - ry[p+q]);
                        /* compute exponential kernel function */
                        rtmp[p] = (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*norm);
                }
                /* insert the diagonal into the list */
                SET_VECTOR_ELT(ans, *rn + q, tmp);
                UNPROTECT(1);
        }
        UNPROTECT(1);

        return(ans);
}

SEXP exponential_kernel_2d_sparse(SEXP x, SEXP y, SEXP l, SEXP var, SEXP n, SEXP m)
{
        R_len_t p, q;
        R_len_t nx;
        R_len_t ny;
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rvar = REAL(var);
        double *rn   = REAL(n);
        double *rm   = REAL(m);
        double *rans, norm;
        double *rtmp;
        SEXP ans, dim, tmp;

        /* check input */
        dim = getAttrib(x, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("x has invalid dimension");
        }
        nx = INTEGER(dim)[0];

        dim = getAttrib(y, R_DimSymbol);
        if(length(dim) != 2 && INTEGER(dim)[1] != 2) {
                error("y has invalid dimension");
        }
        ny = INTEGER(dim)[0];

        if (length(l) != 1) {
                error("l is not a scalar");
        }
        if (length(var) != 1) {
                error("var is not a scalar");
        }
        if (length(n) != 1) {
                error("n is not a scalar");
        }
        if (length(m) != 1) {
                error("m is not a scalar");
        }

        /* define a list that contains the diagonals */
        PROTECT(ans = allocVector(VECSXP, *rn + *rm + 1));
        /* loop over sub-diagonals */
        for(p = *rn; p > 0; p--) {
                /* define a vector that contains a diagonal */
                PROTECT(tmp = allocVector(REALSXP, fmin(nx-p, ny)));
                rtmp = REAL(tmp);
                /* loop through the diagonal */
                for(q = 0; q < fmin(nx-p, ny); q++) {
                        /* compute norm */
                        norm = (rx[p+q + nx*0] - ry[q + ny*0])*(rx[p+q + nx*0] - ry[q + ny*0])
                             + (rx[p+q + nx*1] - ry[q + ny*1])*(rx[p+q + nx*1] - ry[q + ny*1]);
                        /* compute exponential kernel function */
                        rtmp[q] = (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*norm);
                }
                /* insert the diagonal into the list */
                SET_VECTOR_ELT(ans, *rn - p, tmp);
                UNPROTECT(1);
        }
        /* loop over main diagonal and super-diagonals */
        for(q = 0; q <= *rm; q++) {
                /* define a vector that contains a diagonal */
                PROTECT(tmp = allocVector(REALSXP, fmin(nx, ny-q)));
                rtmp = REAL(tmp);
                /* loop through the diagonal */
                for(p = 0; p < fmin(nx, ny-q); p++) {
                        /* compute norm */
                        norm = (rx[p + nx*0] - ry[p+q + ny*0])*(rx[p + nx*0] - ry[p+q + ny*0])
                             + (rx[p + nx*1] - ry[p+q + ny*1])*(rx[p + nx*1] - ry[p+q + ny*1]);
                        /* compute exponential kernel function */
                        rtmp[p] = (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*norm);
                }
                /* insert the diagonal into the list */
                SET_VECTOR_ELT(ans, *rn + q, tmp);
                UNPROTECT(1);
        }
        UNPROTECT(1);

        return(ans);
}
