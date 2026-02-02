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
    fitsfile *afptr, *bfptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, bnaxis, ii,jj,kk, nbin,check = 1, bitpix=-32;
    long npixelsa = 1, npixelsb = 1, firstpixa[2] = {1,1}, firstpixb[2] = {1,1};
    long anaxes[2] = {1,1}, bnaxes[2]={1,1};
    double *apix, *bpix, *opix;
    double biaslevel1=0.,biaslevel2=0.,biaslevel3=0.,biaslevel4=0.,exptime=1.;
    double mjd=0.,rm=RAND_MAX;

    if (argc < 8) { 
      printf("Usage: flatreduce datafile biasfile outfile biaslevel1 biaslevel2 biaslevel3 biaslevel4\n");
      printf("\n");
      printf("Just reduce a flat image\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */
    fits_open_file(&bfptr, argv[2], READONLY, &status);

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_dim(bfptr, &bnaxis, &status);
    fits_get_img_size(afptr, 2, anaxes, &status);
    fits_get_img_size(bfptr, 2, bnaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }

    /* create the new empty output file if the above checks are OK */
    if (check && !fits_create_file(&outfptr, argv[3], &status) )
    {
      /* copy all the header keywords from first image to new output file */
      fits_copy_header(bfptr, outfptr, &status);
      fits_update_key(outfptr, TINT, "BITPIX", &bitpix, NULL, &status);
      fits_delete_key(outfptr, "BZERO", &status);
      if (status!=0) status=0;
      fits_delete_key(outfptr, "BSCALE", &status);
      if (status!=0) status=0;

      fits_read_key(afptr, TFLOAT, "EXPTIME", &exptime, NULL, &status);
      if (status!=0) {
          status=0;
          exptime=1.;
      }
      fits_read_key(afptr, TDOUBLE, "MJD", &mjd, NULL, &status);
      if (status!=0) {
          status=0;
          mjd=0.;
      }
      srand((unsigned long)(mjd*1.e6));

      npixelsa = anaxes[0];  /* no. of pixels to read in each row */
      npixelsb = bnaxes[0];  /* no. of pixels to read in each row */
      nbin = npixelsa/npixelsb;
      /* should have npixelsa >= npixelsb */

      apix = (double *) malloc(npixelsa * sizeof(double)); /* mem for 1 row */
      bpix = (double *) malloc(npixelsb * sizeof(double));
      opix = (double *) malloc(npixelsb * sizeof(double));

      if (apix == NULL || bpix == NULL || opix == NULL) {
        printf("Memory allocation error\n");
        return(1);
      }
      if (argc>5) biaslevel1=atof(argv[4]);
      if (argc>6) biaslevel2=atof(argv[5]);
      if (argc>7) biaslevel3=atof(argv[6]);
      if (argc>8) biaslevel4=atof(argv[7]);

      /* loop over all rows of the plane */
      for (firstpixa[1]=firstpixb[1]=1; firstpixb[1] <= bnaxes[1]; firstpixb[1]++) {
        /* Read both images as doubles, regardless of actual datatype.  */
        /* Give starting pixel coordinate and no. of pixels to read.    */
        /* This version does not support undefined pixels in the image. */

        if ( fits_read_pix(bfptr, TDOUBLE, firstpixb, npixelsb,  NULL, bpix, NULL, &status)  ) break;   /* jump out of loop on error */

        for(ii=0; ii<npixelsb; ii++) opix[ii]=-bpix[ii];
        for (kk=0;kk<nbin;kk++) {
            if ( fits_read_pix(afptr, TDOUBLE, firstpixa, npixelsa, NULL, apix, NULL, &status)) break;

            if (firstpixb[1]<=bnaxes[1]/2) {
               for (jj=0;jj<nbin;jj++) {
                   for(ii=0; ii< npixelsb/2; ii++) opix[ii] += rand()/rm + apix[ii*nbin+jj] - biaslevel1;
                   for(ii=npixelsb/2; ii< npixelsb; ii++) opix[ii] += rand()/rm + apix[ii*nbin+jj] - biaslevel3;
               }
            } else {
               for (jj=0;jj<nbin;jj++) {
                   for(ii=0; ii< npixelsb/2; ii++) opix[ii] += rand()/rm + apix[ii*nbin+jj] - biaslevel2;
                   for(ii=npixelsb/2; ii< npixelsb; ii++) opix[ii] += rand()/rm + apix[ii*nbin+jj] - biaslevel4;
               }
            }

            /*for(ii=0; ii< npixelsb; ii++) {
                for (jj=0;jj<nbin;jj++) opix[ii] += rand()/rm + apix[ii*nbin+jj] - biaslevel;
            }*/
            firstpixa[1]+=1;
        }

        fits_write_pix(outfptr, TDOUBLE, firstpixb, npixelsb, opix, &status); /* write new values to output image */
      }

      fits_close_file(outfptr, &status);
      free(apix);
      free(bpix);
      free(opix);
    }

    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);
 
    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
