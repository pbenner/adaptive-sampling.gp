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
        SEXP ans;
     
        PROTECT(ans = allocMatrix(REALSXP, nx, ny));
        rans = REAL(ans);
        for(i = 0; i < nx; i++) {
                for(j = 0; j < ny; j++) {
                        rans[i + nx*j] = (*rvar)*exp(1.0/(2.0*pow(*rl, 2))*(rx[i] - ry[j])*(rx[i] - ry[j]));
                }
        }
        UNPROTECT(1);

        return(ans);
}
