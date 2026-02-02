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

int main(int argc, char *argv[])
{
    fitsfile *afptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, ii, check = 1, bitpix=-32;
    long npixels = 1, firstpix[2] = {1,1};
    long anaxes[2] = {1,1};
    double *apix, *opix;
    double biaslevel1=0.,biaslevel2=0,biaslevel3=0,biaslevel4=0;
    double mjd=0.,rm=RAND_MAX;

    if (argc < 7) { 
      printf("Usage: darkreduce datafile outfile biaslevel1 biaslevel2 biaslevel3 biaslevel4\n");
      printf("\n");
      printf("Just reduce a dark image\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_size(afptr, 2, anaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }

    /* create the new empty output file if the above checks are OK */
    if (check && !fits_create_file(&outfptr, argv[2], &status) )
    {
      /* copy all the header keywords from first image to new output file */
      fits_copy_header(afptr, outfptr, &status);
      fits_update_key(outfptr, TINT, "BITPIX", &bitpix, NULL, &status);
      fits_delete_key(outfptr, "BZERO", &status);
      if (status!=0) status=0;
      fits_delete_key(outfptr, "BSCALE", &status);
      if (status!=0) status=0;

      fits_read_key(afptr, TDOUBLE, "MJD", &mjd, NULL, &status);
      if (status!=0) {
          status=0;
          mjd=0.;
      }
      srand((unsigned long)(mjd*1.e6));

      npixels = anaxes[0];  /* no. of pixels to read in each row */

      apix = (double *) malloc(npixels * sizeof(double)); /* mem for 1 row */
      opix = (double *) malloc(npixels * sizeof(double));

      if (apix == NULL || opix == NULL) {
        printf("Memory allocation error\n");
        return(1);
      }
      if (argc>3) biaslevel1=atof(argv[3]);
      if (argc>4) biaslevel2=atof(argv[4]);
      if (argc>5) biaslevel3=atof(argv[5]);
      if (argc>6) biaslevel4=atof(argv[6]);

      /* loop over all rows of the plane */
      for (firstpix[1]=1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
        if ( fits_read_pix(afptr, TDOUBLE, firstpix, npixels, NULL, apix, NULL, &status)) break;

        if (firstpix[1]<=anaxes[1]/2) {
               for(ii=0; ii< npixels/2; ii++) opix[ii] = rand()/rm + apix[ii] - biaslevel1;
               for(ii=npixels/2; ii< npixels; ii++) opix[ii] = rand()/rm + apix[ii] - biaslevel3;
        } else {
               for(ii=0; ii< npixels/2; ii++) opix[ii] = rand()/rm + apix[ii] - biaslevel2;
               for(ii=npixels/2; ii< npixels; ii++) opix[ii] = rand()/rm + apix[ii] - biaslevel4;
        }

        /*for (ii=0; ii< npixels; ii++) opix[ii] = rand()/rm + apix[ii] - biaslevel;*/
        fits_write_pix(outfptr, TDOUBLE, firstpix, npixels, opix, &status); /* write new values to output image */
      }

      fits_close_file(outfptr, &status);
      free(apix);
      free(opix);
    }

    fits_close_file(afptr, &status);
 
    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
