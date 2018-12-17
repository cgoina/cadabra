/*
Copyright (C) 2008 Heiko Dankert, California Institute of Technology

 This file is part of QTRAK.

 This program is free software: you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation, either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License
 along with this program.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "mex.h"
#include "matrix.h"

#define NOBIN -1 /* Bin not found */ 

/* Usage: fhistc( *data, *bins, *result ) */

void mexFunction(int nlhs, mxArray *plhs[], int	nrhs, const mxArray *prhs[] )
{
    /* This function takes either single or double data vectors     */
    /* Depending on the type, the appropriate code will be executed */

    register int i, n, k0, k1, bin, nbins;
    register double d;
    register double *dhist = (double *)mxGetPr(prhs[2]);
    register double *ddata = (double *)mxGetData(prhs[0]);
    register double *dedges = (double *)mxGetData(prhs[1]);
    register float s;
    register double *shist = (double *)mxGetPr(prhs[2]);
    register float *sdata = (float *)mxGetData(prhs[0]);
    register float *sedges = (float *)mxGetData(prhs[1]);

    /* Determine the number of data elements and the number of bins */
    /* Then, branch to the appropriate histograming code            */
    
    n = mxGetNumberOfElements(prhs[0]);
    nbins = mxGetNumberOfElements(prhs[1]);
    if( mxIsSingle(prhs[0]) ) goto singlehist;

    /* This is the branch for handling double precision data        */
    
    for( i=0; i < n; i++, ddata++ )
    {
        d = *ddata; 
        bin = NOBIN;

        /* Use a binary search */

        k0 = 0; 
        k1 = nbins - 1;
        if( d >= dedges[0] && d < dedges[nbins-1] )
        {
            bin = (k0+k1)/2;
            while( k0 < k1-1 )
            {
                if ( d >= dedges[bin] ) k0 = bin;
                else k1 = bin;
                bin = (k0+k1)/2;
            }
            bin = k0;
        }

        /* Check for special case */

        if( d == dedges[nbins-1] )
            bin = nbins-1;

        if( bin != NOBIN ) dhist[bin]++;
    }
    
    return;

singlehist:

    /* This is the branch for handling single precision data        */

    for( i=0; i < n; i++, sdata++ )
    {
        s = *sdata; 
        bin = NOBIN;

        /* Use a binary search */
        k0 = 0; 
        k1 = nbins - 1;
        if( s >= sedges[0] && s < sedges[nbins-1] )
        {
            bin = (k0+k1)/2;
            while( k0 < k1-1 )
            {
                if (s >= sedges[bin]) k0 = bin;
                else k1 = bin;
                bin = (k0+k1)/2;
            }
            bin = k0;
        }

        /* Check for special case */
        if( s == sedges[nbins-1] )
            bin = nbins-1;

        if( bin != NOBIN ) shist[bin]++;
    }
}
