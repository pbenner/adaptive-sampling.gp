/* Copyright (C) 2013 Philipp Benner
 *
 * This program is free software; you can redistribute it and/or
 * modify
 * it under the terms of the GNU General Public License as published
 * by
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

SEXP exponential_kernel_1d(SEXP x, SEXP y, SEXP l, SEXP var)
{
        R_len_t i, j;
        R_len_t nx = length(x);
        R_len_t ny = length(y);
        double tmp;
        double *rx   = REAL(x);
        double *ry   = REAL(y);
        double *rl   = REAL(l);
        double *rvar = REAL(var);
        double *rans;
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
                        rans[i + nx*j] =
                                (*rvar)*exp(-1.0/(2.0*(*rl)*(*rl))*(rx[i] - ry[j])*(rx[i] - ry[j]));
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
        double tmp;
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
