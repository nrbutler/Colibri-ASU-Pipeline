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

#include <stdio.h>
#include "fitsio.h"
#include "math.h"

int main(int argc, char *argv[])
{
    fitsfile *afptr, *bfptr, *cfptr, *dfptr, *efptr;  /* FITS file pointers */
    int status = 0, nsat=5000;  /* CFITSIO status value MUST be initialized to zero! */
    int i,k,hough_size=1024;
    long hmax[3] = {0,0,0};
    long imax[3] = {0,0,0};
    long jmax[3] = {0,0,0};
    int m2=hough_size/2,*hpix;
    int anaxis, check = 1, ii;
    long j,nmax=0,npixels = 1, firstpix[2] = {1,1};
    long anaxes[2] = {1,1}, bnaxes[2]={1,1}, cnaxes[2]={1,1}, dnaxes[2]={1,1}, enaxes[2]={1,1};
    float *apix, *bpix, *cpix, *dpix, *epix;
    float tmpf,srcfac=10.,pi=3.14159,wpr,wpi,dth,*cth,*sth,rmax,x,y,fscale,fscale0=1.;
    double nsig=3.,nsig2,v00,v11,diff,vdiff;
    float sat,sat0,gain,gain0,var_min=1.,var_min0=1.;

    if (argc < 6) {
      printf("Usage: weight_clip datafile weightfile maskfile stack_file stack_weight [flxscale] [nsig] [srcfac]\n");
      printf("\n");
      printf("clip out bogux pixels, updates datafile and weightfile\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input images */
    fits_open_file(&bfptr, argv[2], READWRITE, &status);
    fits_open_file(&cfptr, argv[3], READONLY, &status);
    fits_open_file(&dfptr, argv[4], READONLY, &status);
    fits_open_file(&efptr, argv[5], READONLY, &status);

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_size(afptr, 2, anaxes, &status);
    fits_get_img_size(bfptr, 2, bnaxes, &status);
    fits_get_img_size(cfptr, 2, cnaxes, &status);
    fits_get_img_size(dfptr, 2, dnaxes, &status);
    fits_get_img_size(efptr, 2, enaxes, &status);

    if (status) {
       fits_report_error(stderr, status); /* print error message */
       return(status);
    }

    if (anaxis > 3) {
       printf("Error: images with > 3 dimensions are not supported\n");
       check = 0;
    }
         /* check that the input 2 images have the same size */
    else if ( anaxes[0] != bnaxes[0] ||
              anaxes[0] != cnaxes[0] ||
              anaxes[0] != dnaxes[0] ||
              anaxes[0] != enaxes[0] ||
              anaxes[1] != bnaxes[1] ||
              anaxes[1] != cnaxes[1] ||
              anaxes[1] != dnaxes[1] ||
              anaxes[1] != enaxes[1] ) {
       printf("Error: input images don't have same size\n");
       check = 0;
    }

    rmax = 0.5*sqrt(anaxes[0]*anaxes[0]+anaxes[1]*anaxes[1]);

    /* if the above checks are OK */
    if (check==1) {
      /* copy all the header keywords from first image to new output file */
      fits_read_key(afptr, TFLOAT, "FSCALE0", &fscale, NULL, &status);
      if (status!=0) {
          status=0;
          fscale=1.0;
      }
      fits_read_key(afptr, TFLOAT, "SATURATE", &sat, NULL, &status);
      if (status!=0) {
          status=0;
          sat=7000.;
      }

      fits_read_key(dfptr, TFLOAT, "SATURATE", &sat0, NULL, &status);
      if (status!=0) {
          status=0;
          sat0=7000.;
      }

      fits_read_key(afptr, TFLOAT, "GAIN", &gain, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading gain\n");
          gain=1.;
      }
      fits_read_key(afptr, TFLOAT, "VAR0", &var_min, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading VAR0\n");
          var_min=1.;
      }
      fits_read_key(dfptr, TFLOAT, "GAIN", &gain0, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading gain0\n");
          gain0=1.;
      }
      fits_read_key(dfptr, TFLOAT, "VAR0", &var_min0, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading var_min0\n");
          var_min0=1.;
      }

      npixels = anaxes[0];  /* no. of pixels to read in each row */

      apix = (float *) malloc(npixels * sizeof(float)); /* mem for 1 row */
      bpix = (float *) malloc(npixels * sizeof(float));
      cpix = (float *) malloc(npixels * sizeof(float));
      dpix = (float *) malloc(npixels * sizeof(float));
      epix = (float *) malloc(npixels * sizeof(float));
      hpix = (int *) malloc(hough_size*hough_size * sizeof(int));
      sth = (float *) malloc(hough_size * sizeof(float));
      cth = (float *) malloc(hough_size * sizeof(float));

      if (apix == NULL || bpix == NULL || cpix == NULL || dpix == NULL ||
          epix == NULL || hpix == NULL || sth == NULL || cth == NULL ) {
        printf("Memory allocation error\n");
        return(1);
      }
      if (argc>6) fscale0=atof(argv[6]);
      if (argc>7) nsig=atof(argv[7]);
      if (argc>8) srcfac=atof(argv[8]);
      if (fscale0>nsig*fscale) fscale /= fscale0;
      else fscale = 1.0;
      /*fprintf( stderr," nsig=%f srcfac=%f fscale=%f\n",nsig,srcfac,fscale);*/

      for (j=0;j<hough_size*hough_size;j++) hpix[j]=0;

      /*create a grid of sines and cosines between theta=-pi/2 and theta=pi/2*/
      dth=pi/(hough_size-1);
      wpi=sin(dth);
      wpr=sin(0.5*dth); wpr = -2.*wpr*wpr;

      sth[0]=-1; cth[0]=0;
      for (i=1;i<hough_size;i++) {
          sth[i] = (wpr*(tmpf=sth[i-1]) + wpi*cth[i-1]) + sth[i-1];
          cth[i] = (wpr*cth[i-1] - wpi*tmpf) + cth[i-1];
      }

      /*scale gain by exposure map*/
      gain *= var_min;
      gain0 *= var_min0;

      /* loop over all rows of the plane */
      for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
        /* Read both images as floats, regardless of actual datatype.  */
        /* Give starting pixel coordinate and no. of pixels to read.    */
        /* This version does not support undefined pixels in the image. */

        if (fits_read_pix(afptr, TFLOAT, firstpix, npixels, NULL, apix,
                          NULL, &status)  ||
            fits_read_pix(bfptr, TFLOAT, firstpix, npixels,  NULL, bpix,
                          NULL, &status)  ||
            fits_read_pix(cfptr, TFLOAT, firstpix, npixels,  NULL, cpix,
                          NULL, &status)  ||
            fits_read_pix(dfptr, TFLOAT, firstpix, npixels,  NULL, dpix,
                          NULL, &status)  ||
            fits_read_pix(efptr, TFLOAT, firstpix, npixels,  NULL, epix,
                          NULL, &status)  )
            break;   /* jump out of loop on error */

        /*need to work out x and y*/
        y = m2*(firstpix[1]-1-anaxes[1]/2)/rmax;

        for(ii=0; ii< npixels; ii++) {
          if (bpix[ii]>0 && epix[ii]>0 && dpix[ii]<sat0 && apix[ii]<sat) {
              v00 = 1/bpix[ii];
              v11 = 1/epix[ii];
              vdiff = v00 + v11;
              diff = 0.;
              if (apix[ii]>0) {
                  vdiff += v00*apix[ii]/gain;
                  diff = apix[ii];
              }
              if (dpix[ii]>0) {
                  vdiff += dpix[ii]*(dpix[ii]+v11/gain0);
                  diff -= dpix[ii];
              }
              if (diff>nsig*sqrt(vdiff)) {
                  if (ii>20 && firstpix[1]>21 && ii<npixels-20 && firstpix[1]<anaxes[1]-19) {
                      nmax+=1;
                      x = m2*(ii-npixels/2)/rmax;
                      for (i=0;i<hough_size;i++) {
                          j = m2+(long)round(x*cth[i]+y*sth[i]);
                          hpix[j+i*hough_size]+=1;
                      }
                  }
              }
              nsig2=nsig*nsig;
              if (cpix[ii]==0) {
                  nsig2*=srcfac;
                  diff = apix[ii]*fscale-dpix[ii];
              }
              else diff = apix[ii]-dpix[ii];
              if (diff*diff>nsig2*vdiff) bpix[ii]=0;
          }
        }

        fits_write_pix(bfptr, TFLOAT, firstpix, npixels, bpix, &status); /* write new values to output image */
      }

      /* we will consider nsat votes minimum for a strong satellite trail */
      for (k=0;k<3;k++) {
          for (i=0;i<hough_size;i++) {
              for (j=0;j<hough_size;j++) {
                  if (hpix[j+i*hough_size]>hmax[k]) {
                      hmax[k]=hpix[j+i*hough_size];
                      imax[k]=i; jmax[k]=j;
                  }
              }
          }
          /*fprintf( stderr,"Possible sat trail for %s: %ld (nmax=%ld)\n",argv[1],hmax[k],nmax);*/
          if (hmax[k]>nsat) {
               fprintf( stderr,"Removing satellite trail for %s at %ld %ld (votes=%ld/%ld)\n",argv[1],imax[k],jmax[k],hmax[k],nmax);
               for (i=imax[k]-100;i<imax[k]+100;i++) {
                  for (j=jmax[k]-100;j<jmax[k]+100;j++) {
                      if (i>=0 && j>=0 && i<hough_size && j<hough_size) hpix[j+i*hough_size]=0;
                  }
               }
          } else break;
      }

      /* pass over the fits and weight file another time to zero out any satellite trails*/
      if (hmax[0]>nsat) {
          for (firstpix[1] = 1; firstpix[1] <= anaxes[1]; firstpix[1]++) {
                  if (fits_read_pix(afptr, TFLOAT, firstpix, npixels, NULL, apix, NULL, &status) ||
                      fits_read_pix(bfptr, TFLOAT, firstpix, npixels, NULL, bpix, NULL, &status) ) break;
                  y = m2*(firstpix[1]-1-anaxes[1]/2)/rmax;
                  for(ii=0; ii< npixels; ii++) {
                      x = m2*(ii-npixels/2)/rmax;
                      if (hmax[0]>nsat)
                          if ( fabs(m2+x*cth[imax[0]]+y*sth[imax[0]]-jmax[0])<20*m2/rmax ) bpix[ii]=0;
                      if (hmax[1]>nsat)
                          if ( fabs(m2+x*cth[imax[1]]+y*sth[imax[1]]-jmax[1])<20*m2/rmax ) bpix[ii]=0;
                      if (hmax[2]>nsat)
                          if ( fabs(m2+x*cth[imax[2]]+y*sth[imax[2]]-jmax[2])<20*m2/rmax ) bpix[ii]=0;
                  }
                  fits_write_pix(bfptr, TFLOAT, firstpix, npixels, bpix, &status); /* write new values to output image */
          }
      }

      free(apix);
      free(bpix);
      free(cpix);
      free(dpix);
      free(epix);
      free(hpix);
      free(sth);
      free(cth);
    }

    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);
    fits_close_file(cfptr, &status);
    fits_close_file(dfptr, &status);
    fits_close_file(efptr, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
