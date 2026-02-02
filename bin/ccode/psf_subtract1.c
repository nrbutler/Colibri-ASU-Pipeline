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
#include <string.h>
#include <unistd.h>
#include "fitsio.h"
#include "math.h"

void cumsum(float *x, int a, int b) {
    long i,j;

    /* smooth: accumulate vectors for differencing later */
    for (i=1;i<a;i++) x[i]+=x[i-1];
    for (j=a;j<b*a;j+=a) {
      for (i=1;i<a;i++) x[i+j]+=x[i-1+j];
      for (i=0;i<a;i++) x[i+j]+=x[i+j-a];
    }
}


int main(int argc, char *argv[])
{
    fitsfile *afptr, *bfptr, *cfptr, *dfptr, *psffile, *outfptr, *outfptr_wt, *outfptr_s;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int i,j,l,k,a,b,x,y,ix,iy,irad,nbinx,nbiny,srad,rad=-1,filter=1;
    int anaxis, check = 1, sz,nx,ny,irad1,*irad_max,irad_max0;
    long npixels = 1, firstpix[2] = {1,1};
    long i0,j0,j1,j2,j10,j20,anaxes[2] = {1,1}, pnaxes[2] = {1,1};
    float *im, *wim, *ref_im, *ref_wim, *cim, *psf_stack, rad_max, rad0, pmax, wt_min=0.;
    float cim0, vim0, tmp, tmpx, tmpy, dx0, gain, var_min, gain_ref, var_min_ref, sky_lev, sky_lev_ref;
    float *dx, *dy, *numx, *numy, *denomx, *denomy, *denomxy;
    float nmx, nmy, m11, m12, m22;
    float det,delta_x,delta_y,dis2;
    /*float x0,y0;*/

    if (argc < 8) {
      printf("Usage: psf_subtract datafile weightfile reffile weightrefile psfgrid outfile outrmsfile [filter] [smooth radius]\n");
      printf("\n");
      return(0);
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input image */
    fits_open_file(&bfptr, argv[2], READONLY, &status); /* open input weight image */
    fits_open_file(&cfptr, argv[3], READONLY, &status); /* open reference image */
    fits_open_file(&dfptr, argv[4], READONLY, &status); /* open reference weight image */

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
      fits_read_key(cfptr, TFLOAT, "GAIN", &gain_ref, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading gain.\n");
          gain_ref=1.;
      }
      fits_read_key(cfptr, TFLOAT, "VAR0", &var_min_ref, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading VAR0.\n");
          var_min_ref=1.;
      }
      if (status!=0) status=0;
      sky_lev = gain*var_min;
      sky_lev_ref = gain_ref*var_min_ref;
      fprintf(stderr," Sky (gain*var0) levels: %f %f\n",sky_lev,sky_lev_ref);

      a=anaxes[0]; b=anaxes[1];
      npixels = a*b;  /* no. of pixels in image */

      im = (float *) malloc(npixels * sizeof(float)); /* mem for image 1 */
      wim = (float *) malloc(npixels * sizeof(float)); /* mem for image 1 weight */
      ref_im = (float *) malloc(npixels * sizeof(float)); /* mem for image 2 */
      ref_wim = (float *) malloc(npixels * sizeof(float)); /* mem for image 2 weight */
      cim = (float *) malloc(npixels * sizeof(float)); /* mem for output image */

      dx = (float *) malloc(npixels * sizeof(float)); /* mem for intermediate image */
      dy = (float *) malloc(npixels * sizeof(float)); /* mem for intermediate image */
      numx = (float *) malloc(npixels * sizeof(float)); /* mem for intermediate image */
      numy = (float *) malloc(npixels * sizeof(float)); /* mem for intermediate image */
      denomx = (float *) malloc(npixels * sizeof(float)); /* mem for intermediate image */
      denomy = (float *) malloc(npixels * sizeof(float)); /* mem for intermediate image */
      denomxy = (float *) malloc(npixels * sizeof(float)); /* mem for intermediatt image */

      if (im==NULL || wim==NULL || ref_im==NULL || ref_wim==NULL || cim==NULL) {
        printf("Initial memory allocation error\n");
        return(1);
      }
      if (dx==NULL || dy==NULL || numx==NULL || numy==NULL || denomx==NULL || denomy==NULL || denomxy==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      /* now read in the data */
      fits_read_pix(afptr, TFLOAT, firstpix, npixels, NULL, im, NULL, &status);
      fits_read_pix(bfptr, TFLOAT, firstpix, npixels, NULL, wim, NULL, &status);
      fits_read_pix(cfptr, TFLOAT, firstpix, npixels, NULL, ref_im, NULL, &status);
      fits_read_pix(dfptr, TFLOAT, firstpix, npixels, NULL, ref_wim, NULL, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */

      /* get the psf model */
      fits_open_file(&psffile, argv[5], READONLY, &status);
      fits_get_img_size(psffile, 2, pnaxes, &status);
      nx = pnaxes[0]; ny = pnaxes[1];

      fits_read_key(psffile, TINT, "NBINX", &nbinx, NULL, &status);
      fits_read_key(psffile, TINT, "NBINY", &nbiny, NULL, &status);

      sz = nx/nbinx;
      irad = ((sz-1)/2);
      irad_max = (int *) malloc(nbinx*nbiny * sizeof(int));

      psf_stack = (float *) malloc(nx*ny * sizeof(float)); /* mem for psf stack */
      if (psf_stack==NULL) {
        printf("Memory allocation error for psf_stack\n");
        return(1);
      }

      fits_read_pix(psffile, TFLOAT, firstpix, nx*ny, NULL, psf_stack, NULL, &status);

      if (status) fits_report_error(stderr, status); /* print any error message */

      /* prep the images, adding source rms to background rms and inverting */
      for (j0=0;j0<b*a;j0+=a) {
          for (i0=j0;i0<a+j0;i0++) {
              cim[i0]=0;
              if (wim[i0]>wt_min) {
                wim[i0] = 1/wim[i0];
                if (im[i0]>0) wim[i0]*=1+im[i0]/sky_lev;
              } else im[i0]=ref_im[i0]=wim[i0]=ref_wim[i0]=0;
              if (ref_wim[i0]>wt_min) {
                ref_wim[i0] = 1/ref_wim[i0];
                if (ref_im[i0]>0) ref_wim[i0]*=1+ref_im[i0]/sky_lev_ref;
              } else im[i0]=ref_im[i0]=wim[i0]=ref_wim[i0]=0;
          }
      }

      /* book-keeping for psf convolution */

      for (j=0;j<nbiny;j++) {
        for (i=0;i<nbinx;i++) {
          j0 = (i+nx*j)*sz;
          pmax=0;
          for (l=0;l<sz;l++) {
              for (k=0;k<sz;k++) {
                  if (fabs(psf_stack[ k + l*nx + j0 ])>pmax) pmax=fabs(psf_stack[ k + l*nx + j0 ]);
              }
          }
          rad_max=0;
          for (l=0;l<sz;l++) {
              for (k=0;k<sz;k++) {
                  if ( fabs(psf_stack[ k + l*nx + j0 ])>0.01*pmax ) {
                      rad0 = sqrt( (l-irad)*(l-irad)+(k-irad)*(k-irad) );
                      if (rad0>rad_max) rad_max=rad0;
                  }
              }
          }
          irad_max[i+nbinx*j] = (int)ceil(rad_max);
          if (irad_max[i+nbinx*j]>irad) irad_max[i+nbinx*j] = irad;
        }
      }

      /* do the reference image convolution, creating cim from ref_im */
      /* put the total variance in wim */
      for (y=irad;y<b-irad;y++) {
        iy = (y*nbiny)/b; j0 = a*y;
        for (x=irad;x<a-irad;x++) {
          j10 = x+j0;
          if (ref_wim[j10]==0) continue;
          ix = (x*nbinx)/a;
          j20 = ix*sz+irad+nx*(iy*sz+irad);
          cim0=0; vim0=0;
          irad_max0 = irad_max[ix+nbinx*iy];
          for (j=0;j<=2*irad_max0;j++) {
            j1 = j10 + a*(irad_max0-j);
            j2 = j20 + nx*(j-irad_max0);
            irad1 = (int)ceil( sqrt(j*(2*irad_max0-j)) );
            for (i=-irad1;i<=irad1;i++) {
                cim0 += (tmp=psf_stack[i+j2])*ref_im[j1-i];
                vim0 += tmp*tmp*ref_wim[j1-i];
            }
          }
          cim[j10]=cim0; wim[j10]+=vim0; im[j10]-=cim0;
        }
      }

      /* filter difference image for detection purposes (default 1, bad for photometry) */
      if (argc>8) filter = atoi(argv[9]);

      /* this is a smoothing radius, limiting the psf fitting */
      if (argc>9) rad = atoi(argv[10]);
      if (rad<0) rad = irad;

      if (rad>0) {
          fprintf(stderr,"Correcting micro-shifts, smoothing on size %d\n",rad);

          /* calculate psf derivatives for identifying shifts */
          for (j0=0;j0<3*a;j0+=a) for (i0=j0;i0<a+j0;i0++) dx[i0]=dy[i0]=0;
          for (j0=(b-3)*a;j0<b*a;j0+=a) for (i0=j0;i0<a+j0;i0++) dx[i0]=dy[i0]=0;
          for (j0=0;j0<b*a;j0+=a) {
              for (i0=j0;i0<3+j0;i0++) dx[i0]=dy[i0]=0;
              for (i0=a-3+j0;i0<a+j0;i0++) dx[i0]=dy[i0]=0;
          }
          for (j0=3*a;j0<(b-3)*a;j0+=a) {
              for (i0=3+j0;i0<a-3+j0;i0++) {
                dx[i0] = 0.1*( cim[i0-3]-cim[i0+3]-3*(cim[i0-2]-cim[i0+2])+9*(cim[i0-1]-cim[i0+1]) );
              }
              for (i0=3+j0;i0<a-3+j0;i0++) {
                dy[i0] = 0.1*( cim[i0-3*a]-cim[i0+3*a]-3*(cim[i0-2*a]-cim[i0+2*a])+9*(cim[i0-a]-cim[i0+a]) );
              }
          }

          /* calculate the dot-products for fitting shifts */
          for (j0=0;j0<b*a;j0+=a) {
              for (i0=j0;i0<a+j0;i0++) {
                if (wim[i0]>0) {
                    numx[i0] = (tmpx=dx[i0]/wim[i0])*im[i0]; denomx[i0] = tmpx*dx[i0];
                    numy[i0] = (tmpy=dy[i0]/wim[i0])*im[i0]; denomy[i0] = tmpy*dy[i0];
                    denomxy[i0] = tmpy*dx[i0];
                } else numx[i0]=numy[i0]=denomx[i0]=denomy[i0]=denomxy[i0]=0;
              }
          }

          /* smooth: accumulate vectors now, will difference later */
          cumsum(numx,a,b); cumsum(numy,a,b); 
          cumsum(denomx,a,b); cumsum(denomy,a,b); cumsum(denomxy,a,b);

          /* use ref_wim here to track which pixels are fit */
          for (j0=0;j0<a*b;j0+=a) {
              for (i0=j0;i0<a+j0;i0++) ref_wim[i0]=0;
          }

          irad_max0 = rad>irad ? rad : irad;
          for (j0=(irad_max0+1)*a;j0<(b-irad_max0)*a;j0+=a) {
            for (i0=irad_max0+1+j0;i0<a-irad_max0+j0;i0++) {

              if (wim[i0]<=0) continue;
              ref_wim[i0] = 1;

              j1= i0 + irad*a; j2 = i0 - (irad+1)*a;
              /*m11=denomx[irad+j1]-denomx[-irad-1+j1]-denomx[irad+j2]+denomx[-irad-1+j2];
              m22=denomy[irad+j1]-denomy[-irad-1+j1]-denomy[irad+j2]+denomy[-irad-1+j2];*/
              m11=m22=0.;
              for (i=-irad;i<=irad;i++) {
                   j = irad - (int)sqrt(irad*irad-i*i);
                   m11 += denomx[i0+irad-j+a*i]-denomx[i0-irad-1+j+a*i]-denomx[i0+irad-j+a*(i-1)]+denomx[i0-irad-1+j+a*(i-1)];
                   m22 += denomy[i0+irad-j+a*i]-denomy[i0-irad-1+j+a*i]-denomy[i0+irad-j+a*(i-1)]+denomy[i0-irad-1+j+a*(i-1)];
              }
              dx0 = 0.;
              if (m11+m22>0) dx0 += 1.e-2*(m11+m22)/(sz*sz);
              srad = (int)( rad/(1+dx0*wim[i0]) );

              if (srad<1) srad=1;
              if (srad>irad_max0) srad=irad_max0;
              j1= i0 + srad*a; j2 = i0 - (srad+1)*a;

              nmx=numx[srad+j1]-numx[-srad-1+j1]-numx[srad+j2]+numx[-srad-1+j2];
              nmy=numy[srad+j1]-numy[-srad-1+j1]-numy[srad+j2]+numy[-srad-1+j2];
              m11=denomx[srad+j1]-denomx[-srad-1+j1]-denomx[srad+j2]+denomx[-srad-1+j2];
              m22=denomy[srad+j1]-denomy[-srad-1+j1]-denomy[srad+j2]+denomy[-srad-1+j2];
              m12=denomxy[srad+j1]-denomxy[-srad-1+j1]-denomxy[srad+j2]+denomxy[-srad-1+j2];
              m11 += 16.; m22 += 16.; /*priors, 0.25 pixel rms shift in x and y*/

              det = m11*m22-m12*m12;
              if (det>0) {
                delta_x = (m22*nmx-m12*nmy)/det;
                delta_y = (m11*nmy-m12*nmx)/det;
                dis2 = delta_x*delta_x+delta_y*delta_y;
                if (dis2<1) {
                    im[i0] -= delta_x*dx[i0] + delta_y*dy[i0];
                    cim[i0] += delta_x*dx[i0] + delta_y*dy[i0];
                }
              }
              cim0 = fabs(cim[i0]);
              if (cim0>2*sqrt(wim[i0])) dx0 += cim0;
              if (filter==1) im[i0] /= sqrt(1+dx0*dx0);
              /*if (filter==1) {
                  dx0 = (1+2*dx0)/(1+dx0)/(1+dx0);
                  im[i0] *= dx0;
                  wim[i0] /= dx0*dx0;
              }*/
            }
          }

      }

      /* finally, make sure we invert and zero out when weight is 0 */
      for (j0=0;j0<a*b;j0+=a) {
          for (i0=j0;i0<a+j0;i0++) {
              if (wim[i0]>0 && ref_wim[i0]>0) {
                  wim[i0]=1./wim[i0];
                  if (im[i0]>0) wim[i0]*=1+im[i0]/sky_lev;
              }
              else wim[i0]=im[i0]=0;
          }
      }

      /* write out the results */
      fits_create_file(&outfptr, argv[6], &status);
      fits_create_file(&outfptr_wt, argv[7], &status);
      fits_create_file(&outfptr_s, argv[8], &status);
 
      fits_copy_header(afptr, outfptr, &status);
      fits_copy_header(afptr, outfptr_wt, &status);
      fits_copy_header(afptr, outfptr_s, &status);

      fits_write_pix(outfptr, TFLOAT, firstpix, npixels, im, &status);
      fits_write_pix(outfptr_wt, TFLOAT, firstpix, npixels, wim, &status);
      fits_write_pix(outfptr_s, TFLOAT, firstpix, npixels, cim, &status);
 
      fits_close_file(outfptr, &status);
      fits_close_file(outfptr_wt, &status);
      fits_close_file(outfptr_s, &status);

      free(im);
      free(wim);
      free(ref_im);
      free(ref_wim);
      free(cim);

      free(dx);
      free(dy);
      free(numx);
      free(numy);
      free(denomx);
      free(denomy);
      free(denomxy);

      free(psf_stack);
      free(irad_max);
    }
    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);
    fits_close_file(cfptr, &status);
    fits_close_file(dfptr, &status);
    fits_close_file(psffile, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
