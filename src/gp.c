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

SEXP exponential_kernel_1d(SEXP x, SEXP y)
{
        R_len_t i, j, nx = length(x), ny = length(y);
        double tmp, *rx = REAL(x), *ry = REAL(y), *rans;
        SEXP ans, dim, dimnames;
     
        PROTECT(ans = allocVector(REALSXP, nx*ny));
        rans = REAL(ans);
        for(i = 0; i < nx; i++) {
                tmp = rx[i];
                for(j = 0; j < ny; j++)
                        rans[i + nx*j] = tmp * ry[j];
        }
     
        PROTECT(dim = allocVector(INTSXP, 2));
        INTEGER(dim)[0] = nx; INTEGER(dim)[1] = ny;
        setAttrib(ans, R_DimSymbol, dim);
     
        PROTECT(dimnames = allocVector(VECSXP, 2));
        SET_VECTOR_ELT(dimnames, 0, getAttrib(x, R_NamesSymbol));
        SET_VECTOR_ELT(dimnames, 1, getAttrib(y, R_NamesSymbol));
        setAttrib(ans, R_DimNamesSymbol, dimnames);
     
        UNPROTECT(3);
        return(ans);
}
