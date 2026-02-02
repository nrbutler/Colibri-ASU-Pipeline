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

#define MAX(a,b) ((a) > (b) ? (a) : (b))

/* function to calculate a bunch of aperture fluxes */
int calc_aper(long nx, long ny, long nk, float *im, float *wim, long *seg, float *bim, float *mx, float *my, float raper, float gain, float *flx, float *dflx, int *flags, int ap) {

    long k;
    int x,y,x1,y1,mx0,my0,xmin,xmax,ymin,ymax,ignore,flagged;
    float la,dx,dy,dys,r,rmin,rmax,a,ai;
    double flx0,vflx0,pix,pix0,var0;

    for (k=0;k<nk;k++) {
        xmin = (int)(mx[k]-raper-0.500001);
        xmax = (int)(mx[k]+raper+0.499999);
        ymin = (int)(my[k]-raper-0.500001);
        ymax = (int)(my[k]+raper+0.499999);

        if (xmin<1 || xmax>=nx-1 || ymin<1 || ymax>=ny-1) {
            flx[k*2+ap]=0; dflx[k*2+ap]=0;
            continue;
        }
        mx0 = (int)ceil(mx[k]-1);
        my0 = (int)ceil(my[k]-1);
        if (wim[mx0+nx*my0]<=0) {
            flx[k*2+ap]=0; dflx[k*2+ap]=0;
            continue;
        }

        flx0=0; vflx0=0; flagged=0;
        for (y=ymin;y<=ymax;y++) {
            dy = y-my[k]+1;
            dys = dy*dy;
            for (x=xmin;x<=xmax;x++) {
                dx = x-mx[k]+1;

                ai = MAX(fabs(dx),fabs(dy));
                a = (ai>0) ? 0.5/ai : 0;

                r = sqrt(dys+dx*dx);
                rmin = r*(1-a);

                if (rmin>=raper) continue;

                rmax = r*(1+a);
                la = (raper<rmax) ? (raper-rmin)/(rmax-rmin) : 1;

                ignore=1;
                if (seg[x+nx*y]==0 || seg[x+nx*y]==k+1) ignore=0;

                pix=0; var0=0;
                if (wim[x+nx*y]>0 && ignore==0) {
                    pix=im[x+nx*y]-bim[x+nx*y];
                    var0=1/wim[x+nx*y];
                } else {
                    x1 = (long)(2*(mx[k]-1)+0.49999-x);
                    y1 = (long)(2*(my[k]-1)+0.49999-y);
                    ignore=1;
                    if (seg[x1+nx*y1]==k+1) ignore=0;
                    if (wim[x1+nx*y1]>0 && ignore==0) {
                        pix=im[x1+nx*y1]-bim[x1+nx*y1];
                        var0=1/wim[x1+nx*y1];
                    }
                }

                if (var0>0) {
                    flx0 += la*pix;
                    pix0 = 0;
                    if (pix>0) pix0 = pix;
                    vflx0 += la*la*var0*(1+pix0/gain);
                }
                if (wim[x+nx*y]==0) flagged=1;
            }
        }
        flx[k*2+ap] = flx0; dflx[k*2+ap] = sqrt(vflx0);
        flags[k] = flagged;
    }
    return 0;
}

int main(int argc, char *argv[])
{
    fitsfile *afptr, *bfptr, *cfptr, *dfptr, *efptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, check = 1, anynulls;
    long npixels = 1, firstpix[2] = {1,1}, nrows;
    long k,k1,anaxes[2] = {1,1};
    float *apix, *dpix, *bpix, *mx, *my, *flx, *dflx;
    long *cpix;
    int x0,y0,*flags_initial,*flags_final;
    float x,y,raper1,raper2,gain=1.,var_min=1.,phot_diam=6.,fwhm,exptime;
    int colnum_x,colnum_y,colnum_f,colnum_df,colnum_flags;
    char dobs[FLEN_VALUE], dobe[FLEN_VALUE];

    if (argc < 6) {
      printf("Usage: calc_phot datafile weightfile segfile backfile catfile [crpix1] [crpix2]\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input image */
    fits_open_file(&bfptr, argv[2], READONLY, &status); /* open input weight image */
    fits_open_file(&cfptr, argv[3], READONLY, &status); /* open input image segmentation map */
    fits_open_file(&dfptr, argv[4], READWRITE, &status); /* open background map */
    fits_open_file(&efptr, argv[5], READWRITE, &status); /* open catalog table */
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

    if (check) {
      /* need to update some header keywords */
      fits_read_key(afptr, TFLOAT, "FWHM", &fwhm, NULL, &status);
      fits_read_key(afptr, TFLOAT, "EXPTIME", &exptime, NULL, &status);
      fits_read_key(afptr, TSTRING, "DATE-OBS", dobs, NULL, &status);
      fits_read_key(afptr, TSTRING, "DATE-OBE", dobe, NULL, &status);
      if (status!=0) status=0;

      /* get some kewords */
      fits_read_key(afptr, TFLOAT, "GAIN", &gain, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading gain.\n");
          gain=1.;
      }
      fits_read_key(afptr, TFLOAT, "VAR0", &var_min, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading VAR0.\n");
          var_min=1.;
      }
      fits_read_key(afptr, TFLOAT, "PHOTDIAM", &phot_diam, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading phot_diam.\n");
          phot_diam=6.;
      }
      fits_read_key(afptr, TFLOAT, "CRPIX1", &x, NULL, &status);
      fits_read_key(afptr, TFLOAT, "CRPIX2", &y, NULL, &status);
      if (status!=0) status=0;

      dobs[10] = '_'; dobe[10] = '_';
      k1=0;
      for (k=0;k<23;k++) {
          dobs[k-k1] = dobs[k];
          dobe[k-k1] = dobe[k];
          if (k==4 || k==7 || k==13 || k==16) k1+=1;
      }
      for (k=19;k<23;k++) dobs[k]=dobe[k]=' ';

      fits_update_key(efptr, TFLOAT, "FWHM", &fwhm, NULL, &status);
      fits_update_key(efptr, TFLOAT, "EXPTIME", &exptime, NULL, &status);
      fits_update_key(efptr, TSTRING, "T0", dobs, NULL, &status);
      fits_update_key(efptr, TSTRING, "T1", dobe, NULL, &status);
      if (status!=0) status=0;

      raper1 = phot_diam/2.;
      raper2 = 3*raper1;

      npixels = anaxes[0]*anaxes[1];  /* no. of pixels in image */

      x0=0; y0=0;
      if (argc>6) x0=(int)(atof(argv[6])-x);
      if (argc>7) y0=(int)(atof(argv[7])-y);

      apix = (float *) malloc(npixels * sizeof(float)); /* mem for image 1 */
      bpix = (float *) malloc(npixels * sizeof(float)); /* mem for image 1 */
      cpix = (long *) malloc(npixels * sizeof(long)); /* mem for image 2 */
      dpix = (float *) malloc(npixels * sizeof(float)); /* mem for image 2 */

      if (apix == NULL || bpix==NULL || cpix==NULL || dpix==NULL ) {
        printf("Memory allocation error\n");
        return(1);
      }
      fits_read_pix(afptr, TFLOAT, firstpix, npixels, NULL, apix, NULL, &status);
      fits_read_pix(bfptr, TFLOAT, firstpix, npixels, NULL, bpix, NULL, &status);
      fits_read_pix(cfptr, TLONG, firstpix, npixels, NULL, cpix, NULL, &status);
      fits_read_pix(dfptr, TFLOAT, firstpix, npixels, NULL, dpix, NULL, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */

      /* get table size */
      fits_movrel_hdu(efptr, 1, NULL, &status);
      fits_get_num_rows(efptr, &nrows, &status);

      /* find the columns numbers */
      fits_get_colnum(efptr, CASEINSEN, "x_image", &colnum_x, &status);
      fits_get_colnum(efptr, CASEINSEN, "y_image", &colnum_y, &status);
      fits_get_colnum(efptr, CASEINSEN, "flux_aper", &colnum_f, &status);
      fits_get_colnum(efptr, CASEINSEN, "fluxerr_aper", &colnum_df, &status);
      fits_get_colnum(efptr, CASEINSEN, "flags", &colnum_flags, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */

      mx = (float *) malloc(nrows * sizeof(float));
      my = (float *) malloc(nrows * sizeof(float));
      flx = (float *) malloc(2*nrows * sizeof(float));
      dflx = (float *) malloc(2*nrows * sizeof(float));
      flags_initial = (int *) malloc(nrows * sizeof(int));
      flags_final = (int *) malloc(nrows * sizeof(int));
      if (mx == NULL || my==NULL || flx==NULL || dflx==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      /* read in the table data */
      fits_read_col(efptr, TFLOAT, colnum_x, 1, 1, nrows, NULL, mx, &anynulls, &status);
      fits_read_col(efptr, TFLOAT, colnum_y, 1, 1, nrows, NULL, my, &anynulls, &status);
      fits_read_col(efptr, TINT, colnum_flags, 1, 1, nrows, NULL, flags_initial, &anynulls, &status);
      /* these are length 2 vectors */
      fits_read_col(efptr, TFLOAT, colnum_f, 1, 1, nrows*2, NULL, flx, &anynulls, &status);
      fits_read_col(efptr, TFLOAT, colnum_df, 1, 1, nrows*2, NULL, dflx, &anynulls, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */

      /* shift the mx,my if necessary */
      if (x0!=0 || y0!=0) {
          for (k=0;k<nrows;k++) {
              mx[k] -= x0; my[k] -= y0;
          }
      }

      /* get the aperture fluxes from the images */
      fprintf(stderr," Extracting flux for %s in %.1f pixel radius aperture...\n",argv[1],raper1);
      calc_aper(anaxes[0], anaxes[1], nrows, apix, bpix, cpix, dpix, mx, my, raper1, gain*var_min, flx, dflx, flags_final, 0);

      fprintf(stderr," Extracting flux for %s in %.1f pixel radius aperture...\n",argv[1],raper2);
      calc_aper(anaxes[0], anaxes[1], nrows, apix, bpix, cpix, dpix, mx, my, raper2, gain*var_min, flx, dflx, flags_final, 1);

      /* use flags_initial for deblending and flags_final for edges */
      for (k=0;k<nrows;k++) {
          if (flags_initial[k]>0 && flags_final[k]==0) flags_final[k] = 2;
      }

      /* now write out the new columns */
      /*fits_write_col(efptr, TINT, colnum_flags, 1, 1, nrows, flags_final, &status);*/
      fits_write_col(efptr, TFLOAT, colnum_f, 1, 1, nrows*2, flx, &status);
      fits_write_col(efptr, TFLOAT, colnum_df, 1, 1, nrows*2, dflx, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */

      free(apix);
      free(bpix);
      free(cpix);
      free(dpix);
      free(mx);
      free(my);
      free(flx);
      free(dflx);
      free(flags_initial);
      free(flags_final);
    }
    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);
    fits_close_file(cfptr, &status);
    fits_close_file(dfptr, &status);
    fits_close_file(efptr, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
