/*
  This file was downloaded from the CFITSIO utilities web page:
    http://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html

  That page contains this text:
    You may freely modify, reuse, and redistribute these programs as you wish.

  We assume it was originally written by the CFITSIO authors (primarily William
  D. Pence).

  We (the Astrometry.net team) have modified it slightly.
  # Licensed under a 3-clause BSD style license - see LICENSE


*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include "fitsio.h"
/* note #undef's at end of file */
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

double select_med(unsigned long k, unsigned long n, double arr[])
{
	unsigned long i,ir,j,l,mid;
	double a,temp;

	l=1;
	ir=n;
	for (;;) {
		if (ir <= l+1) {
			if (ir == l+1 && arr[ir] < arr[l]) {
				SWAP(arr[l],arr[ir])
			}
			return arr[k];
		} else {
			mid=(l+ir) >> 1;
			SWAP(arr[mid],arr[l+1])
			if (arr[l] > arr[ir]) {
				SWAP(arr[l],arr[ir])
			}
			if (arr[l+1] > arr[ir]) {
				SWAP(arr[l+1],arr[ir])
			}
			if (arr[l] > arr[l+1]) {
				SWAP(arr[l],arr[l+1])
			}
			i=l+1;
			j=ir;
			a=arr[l+1];
			for (;;) {
				do i++; while (arr[i] < a);
				do j--; while (arr[j] > a);
				if (j < i) break;
				SWAP(arr[i],arr[j])
			}
			arr[l+1]=arr[j];
			arr[j]=a;
			if (j >= k) ir=j-1;
			if (j <= k) l=i;
		}
	}
}

int main(int argc, char *argv[])
{
    fitsfile *afptr;
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, check = 1, ii;
    long npixels = 1, firstpix[2] = {1,1};
    long anaxes[2] = {1,1};
    double *apix,*allpix,mjd=0.,rm=RAND_MAX;

    if (argc != 2) { 
      printf("Usage: immed datafile\n");
      printf("\n");
      printf("Just give the median of an image\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */
    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_size(afptr, 2, anaxes, &status);

    fits_read_key(afptr, TDOUBLE, "MJD", &mjd, NULL, &status);
    if (status!=0) {
        status=0;
        mjd=0.;
    }
    srand((unsigned long)(mjd*1.e6));

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }

    if (check)
    {
      npixels = anaxes[0];  /* no. of pixels to read in each row */

      apix = (double *) malloc(npixels * sizeof(double)); /* mem for 1 row */
      allpix = (double *) malloc(npixels*anaxes[1] * sizeof(double)); /* mem for 1 row */

      if (apix == NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      /* loop over all rows of the plane */
      for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++)
      {
        /* Read both images as doubles, regardless of actual datatype.  */
        /* Give starting pixel coordinate and no. of pixels to read.    */
        /* This version does not support undefined pixels in the image. */

        if (fits_read_pix(afptr, TDOUBLE, firstpix, npixels, NULL, apix, NULL, &status)) break;   /* jump out of loop on error */
        for(ii=0; ii< npixels; ii++) allpix[ii+(firstpix[1]-1)*anaxes[0]]=rand()/rm+apix[ii];

      }

      free(apix);
    }

    fits_close_file(afptr, &status);

    printf("%.8f\n",select_med(anaxes[0]*anaxes[1]/2, anaxes[0]*anaxes[1], allpix));

    free(allpix);
 
    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}

#undef SWAP
