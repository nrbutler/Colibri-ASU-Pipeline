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

static int iminarg1,iminarg2;
#define IMIN(a,b) (iminarg1=(a),iminarg2=(b),(iminarg1) < (iminarg2) ? (iminarg1) : (iminarg2))
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp;

float median(unsigned long n, float arr[]) {
	unsigned long i,ir,j,l,mid,k=n/2;
	float a,temp;

	l=0;
	ir=n-1;
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

/* function to refine source positions */
int refine_mxy(int nx, int ny, long nk, float *im, long *idx, float *mx, float *my, float *flx, float *dflx, float fwhm) {

    unsigned long k,k0;
    int i,j,dj,i0,j0,mx0,my0;
    float s0[8],s1,delta_x,delta_y,psfx,psfy,dx,dy,maxx,dmx_best,dmy_best,fac;

    /*printf("%d %d %ld\n",nx,ny,nk);*/
    if (fwhm<2) fac = 0.6931;
    else fac = 4*log(2)/(fwhm*fwhm);

    for (k0=0;k0<nk;k0++) {
        k=idx[k0];

        if (flx[2*k]<=dflx[2*k]) continue;

        mx0 = (int)ceil(mx[k]-1); dx = mx0-mx[k]+1;
        my0 = (int)ceil(my[k]-1); dy = my0-my[k]+1;

        maxx=0;
        dmx_best = 0.;
        dmy_best = 0.;
        for (i=0;i<=20;i++) {   
            for (j0=-4;j0<4;j0++) {
                s0[j0+4]=0;
                for (i0=-4;i0<4;i0++) {
                  delta_x = i0+dx-0.1*i+1.;
                  psfx = exp(-fac*delta_x*delta_x);
                  s0[j0+4] += im[mx0+i0+nx*(my0+j0)]*psfx;
                }
            }
            dj = (int)sqrt(i*(20-i));
            for (j=10-dj;j<=10+dj;j++) {   
              s1=0.;
              for (j0=-4;j0<4;j0++) {
                delta_y = j0+dy-0.1*j+1.;
                psfy = exp(-fac*delta_y*delta_y);
                s1 += s0[j0+4]*psfy;
              }
              if (s1>0 && s1>maxx) {
                dmx_best=0.1*i-1.;
                dmy_best=0.1*j-1.;
                maxx=s1;
              }
            }
        }
        mx[k] += dmx_best;
        my[k] += dmy_best;
    }
    return 0;
}

int calc_ap(int nx, int ny, long nk, float *im, float *wim, long *seg, long *idx, float *mx, float *my, float *flx, float *dflx, float raper0, float gain, int *flags) {
    /* local variables */
    int x,y,x1,y1,ymin,ymax,xmin,xmax,flagged;
    int mx0,my0,ignore;
    unsigned long k,k0;
    float dx,dy,dys,r,rmin,rmax;
    float pix,pix0,a,ai,var0;
    float la,la0,flx0,vflx0,flx1,vflx1;
    float raper=3*raper0;

    /* loop over the psf stars */
    for (k0=0;k0<nk;k0++) {
        k=idx[k0];
        flagged=1;

        xmin = (int)(mx[k]-raper-0.500001);
        xmax = (int)(mx[k]+raper+0.499999);
        ymin = (int)(my[k]-raper-0.500001);
        ymax = (int)(my[k]+raper+0.499999);

        if (xmin<0 || xmax>nx-1 || ymin<0 || ymax>ny-1 || nx<=0 || ny<=0) {
            flx[k*2]=dflx[k*2]=flx[k*2+1]=dflx[k*2+1]=0;
            flags[k] = flagged;
            continue;
        }

        mx0 = (int)ceil(mx[k]-1); dx = mx0-mx[k]+1;
        my0 = (int)ceil(my[k]-1); dy = my0-my[k]+1;
        if (wim[mx0+nx*my0]<=0) {
            flx[k*2]=dflx[k*2]=flx[k*2+1]=dflx[k*2+1]=0;
            flags[k] = flagged;
            continue;
        }

        /* calculate aperture fluxes */
        flx0=0; vflx0=0; flx1=0; vflx1=0; flagged=0;
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

                ignore=1;
                if (seg[x+nx*y]==0 || seg[x+nx*y]==k+1) ignore=0;

                pix=0; var0=0;
                if (wim[x+nx*y]>0 && ignore==0) {
                    pix=im[x+nx*y];
                    var0=1/wim[x+nx*y];
                } else {
                    x1 = (int)(2*(mx[k]-1)+0.49999-x);
                    y1 = (int)(2*(my[k]-1)+0.49999-y);
                    if (x1>=0 && x1<nx && y1>=0 && y1<ny) {
                        ignore=1;
                        if (seg[x1+nx*y1]==k+1) ignore=0;
                        if (wim[x1+nx*y1]>0 && ignore==0) {
                            pix=im[x1+nx*y1];
                            var0=1/wim[x1+nx*y1];
                        }
                    }
                }

                if (var0>0) {
                    rmax = r*(1+a);
                    la = (raper<rmax) ? (raper-rmin)/(rmax-rmin) : 1;
                    la0 = (raper0<rmax) ? (raper0-rmin)/(rmax-rmin) : 1;

                    /* add up the fluxes */
                    if (la>0) {
                        pix0 = pix>0 ? pix : 0;
                        if (la0>0) {
                            flx0 += la0*pix;
                            vflx0 += la0*la0*var0*(1+pix0/gain);
                        }
                        flx1 += la*pix;
                        vflx1 += la*la*var0*(1+pix0/gain);
                    }
                }
                if (wim[x+nx*y]==0) flagged=1;
            }
        }

        flx[2*k] = flx0;
        dflx[2*k] = sqrt(vflx0);
        flx[2*k+1] = flx1;
        dflx[2*k+1] = sqrt(vflx1);
        flags[k] = flagged;
    }

    return 0;
}

/* function to evaulate a psf grid on the data */
int calc_psf(int nx, int ny, long nk, float *im, float *wim, long *seg, float *psf, long *idx, float *mx, float *my, float *flx, float *dflx, float raper0, float gain, float spsf) {
    /* local variables */
    int x,y,x1,y1,ymin,ymax,xmin,xmax,ki,kj,imin,imax,jmin,jmax;
    int xp,yp,sz,mx0,my0,ignore,irad;
    unsigned long k,k0,jj1;
    float dx,dy,dx1,dy1,dys,r,rmin,rmax;
    float pix,pix0,a,ai,tmp,var0,wt,psf0,gain_eff;
    float *L4dx,*L4dy,la0,flx0,vflx0,norm,flux,dflux;

    L4dx = (float *) malloc(8 * sizeof(float));
    L4dy = (float *) malloc(8 * sizeof(float));

    irad = (int)ceil(2*raper0);
    sz = 1+2*irad;

    /* loop over the psf stars */
    for (k0=0;k0<nk;k0++) {
        k=idx[k0];

        xmin = (int)(mx[k]-3*raper0-0.500001);
        xmax = (int)(mx[k]+3*raper0+0.499999);
        ymin = (int)(my[k]-3*raper0-0.500001);
        ymax = (int)(my[k]+3*raper0+0.499999);

        if (xmin<0 || xmax>nx-1 || ymin<0 || ymax>ny-1 || nx<=0 || ny<=0) continue;

        mx0 = (int)ceil(mx[k]-1); dx = mx0-mx[k]+1;
        my0 = (int)ceil(my[k]-1); dy = my0-my[k]+1;
        if (wim[mx0+nx*my0]<=0) continue;

        /* make the psf interpolation function */
        for(ki=0;ki<8;ki++) {
            dx1 = dx-ki+3;
            dy1 = dy-ki+3;
            L4dx[ki] = (dx1==0) ? 1 : sin(3.14159*dx1)*sin(0.7854*dx1)/(2.4674*dx1*dx1);
            L4dy[ki] = (dy1==0) ? 1 : sin(3.14159*dy1)*sin(0.7854*dy1)/(2.4674*dy1*dy1);
        }

        /* calculate psf fluxes */
        flx0=0; vflx0=0; norm=0;
        for (y=ymin;y<=ymax;y++) {
            dy = y-my[k]+1;
            dys = dy*dy;
            yp = y-my0+irad;
            jmin=-IMIN(3,yp); jmax=IMIN(4,2*irad-yp);

            for (x=xmin;x<=xmax;x++) {
                dx = x-mx[k]+1;

                ai = MAX(fabs(dx),fabs(dy));
                a = (ai>0) ? 0.5/ai : 0;

                r = sqrt(dys+dx*dx);
                rmin = r*(1-a);

                if (rmin>=raper0) continue;

                ignore=1;
                if (seg[x+nx*y]==0 || seg[x+nx*y]==k+1) ignore=0;

                pix=0; var0=0;
                if (wim[x+nx*y]>0 && ignore==0) {
                    pix=im[x+nx*y];
                    var0=1/wim[x+nx*y];
                } else {
                    x1 = (int)(2*(mx[k]-1)+0.49999-x);
                    y1 = (int)(2*(my[k]-1)+0.49999-y);
                    if (x1>=0 && x1<nx && y1>=0 && y1<ny) {
                        ignore=1;
                        if (seg[x1+nx*y1]==k+1) ignore=0;
                        if (wim[x1+nx*y1]>0 && ignore==0) {
                            pix=im[x1+nx*y1];
                            var0=1/wim[x1+nx*y1];
                        }
                    }
                }

                if (var0>0) {
                    xp = x-mx0+irad;

                    rmax = r*(1+a);
                    la0 = (raper0<rmax) ? (raper0-rmin)/(rmax-rmin) : 1;

                    /* map the psf to the image grid */
                    if (la0>0) {
                        imin=-IMIN(3,xp); imax=IMIN(4,2*irad-xp);
                        psf0=0;
                        for(kj=jmin;kj<=jmax;kj++) {
                          jj1 = xp + sz*(yp + kj);
                          for(ki=imin;ki<=imax;ki++) psf0 += L4dx[ki+3]*L4dy[kj+3]*psf[jj1+ki];
                        }
                        if (psf0>0) {
                            dflux = flux = 0;
                            if (flx[2*k]>0) {
                                dflux = flx[2*k];
                                flux = dflux*psf0;
                            }
                            if (dflx[2*k]>1) dflux /= dflx[2*k];
                            gain_eff = gain/(1+0.1*dflux);
                            wt = la0/var0/( 1 + flux/gain_eff );
                            tmp = psf0*wt;
                            flx0 += tmp*pix;
                            norm += tmp*psf0;
                            pix0 = pix>0 ? pix : 0;
                            vflx0 += tmp*tmp*var0*(1+pix0/gain);
                        }
                    }
                }
            }
        }
        if (norm>0) {
            flx[2*k] = flx0*spsf/norm;
            dflx[2*k] = sqrt(vflx0)*spsf/norm;
        } else {
            fprintf(stderr,"Warning: negative or zero norm!\n");
        }
    }
    free(L4dx);
    free(L4dy);

    return 0;
}

/* function to calculate a grid of point-spread functions (psfs) */
int build_median_psf(int nx, int ny, long npsf_stars, float *im, float *wim, long *seg, float *psf, long *idx, float *mx, float *my, float *flx, int irad, float gain) {
    /* local variables */
    int k0,k00,i,j,ir,jr,sz,sz2,ki,kj,mx0,my0,norm0;
    long k,j1;
    float dx,dy,dx1,dy1,*L4dx,*L4dy,im0,im1,im2,*psf_ar,*sar,*norm;

    sz = 1+2*irad;
    sz2 = sz*sz;

    L4dx = (float *) malloc(8 * sizeof(float));
    L4dy = (float *) malloc(8 * sizeof(float));
    psf_ar = (float *) malloc(npsf_stars*sz2 * sizeof(float));
    sar = (float *) malloc(npsf_stars * sizeof(float));
    norm = (float *) malloc(npsf_stars * sizeof(float));

    /* loop over the psf stars */
    for (k0=0;k0<npsf_stars;k0++) {
        k = idx[k0];
        norm[k0]=0;

        mx0 = (int)ceil(mx[k]-1); dx = mx0-mx[k]+1;
        my0 = (int)ceil(my[k]-1); dy = my0-my[k]+1;

        for (j=0;j<sz;j++) {
           j1 = j*sz+sz2*k0;
           for (i=0;i<sz;i++) psf_ar[i+j1]=0;
        }

        if ( mx0<irad || mx0+irad>nx-1 || my0<irad || my0+irad>ny-1 ) continue;
 
        /* make the interpolation function */
        for(ki=0;ki<8;ki++) {
            dx1 = dx-ki+3;
            dy1 = dy-ki+3;
            L4dx[ki] = (dx1==0) ? 1 : sin(3.14159*dx1)*sin(0.7854*dx1)/(2.4674*dx1*dx1);
            L4dy[ki] = (dy1==0) ? 1 : sin(3.14159*dy1)*sin(0.7854*dy1)/(2.4674*dy1*dy1);
        }

        /* map the psf star to the correct (centered) grid */
        for (j=-irad;j<=irad;j++) {
            j1 = mx0 + nx*(j + my0);
            jr = j+irad;
            for (i=-irad;i<=irad;i++) {
                if ( seg[i+j1]==k+1 || seg[i+j1]==0 ) {
                    ir = i+irad;
                    if (wim[i+j1]>0) {
                        im0 = im[i+j1];
                        for(kj=-IMIN(3,jr);kj<=IMIN(4,2*irad-jr);kj++) {
                          im1 = im0*L4dy[kj+3];
                          for(ki=-IMIN(3,ir);ki<=IMIN(4,2*irad-ir);ki++) {
                            im2 = im1*L4dx[ki+3];
                            psf_ar[ir+ki+sz*(jr+kj)+k0*sz2] += im2;
                            norm[k0] += im2;
                          }
                        }
                    }
                }
            }
        }

    }

    /* now get the median at each pixel */
    norm0=0;
    for (j=0;j<sz;j++) {
        for (i=0;i<sz;i++) {
            j1 = i+j*sz;
            k00=0;
            for (k0=0;k0<npsf_stars;k0++) {
                if (norm[k0]>0 && flx[2*idx[k0]]>0) {
                    norm0+=1;
                    sar[k00] = psf_ar[j1+k0*sz2]/flx[2*idx[k0]];
                    k00++;
                }
            }
            if (k00>0) psf[j1] = median(k00,sar);
            else psf[j1]=0;
        }
    }

    free(L4dx);
    free(L4dy);
    free(psf_ar);
    free(sar);
    free(norm);

    return norm0;
}


int main(int argc, char *argv[])
{
    fitsfile *afptr, *bfptr, *cfptr, *dfptr, *efptr, *outfptr;  /* FITS file pointers */
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int anaxis, check = 1, anynulls, make_psf, make_calc=1, nx0, ny0, sz,smooth_iter=1, refine_pos=0;
    int npix_x=0,npix_y=0,npsf=15;
    long npixels = 1, firstpix[2] = {1,1}, lastpix[2] = {1,1}, inc[2] = {1,1}, nrows;
    long k,k0,j1,jj1,anaxes[2] = {1,1}, enaxes1[2] = {1,1}, npsf_stars=0;
    float *apix, *bpix, *mx, *my, *flx, *dflx, *psf_stack, *psf;
    long *cpix, *idx, *idx_ar, *idx_ar0, *nixy, *cnixy, *nixy_iter, *nixy0, *cnixy0, *nixy_iter0;
    int x0,y0,irad,nbinx,nbiny,nx,ny,*flags_initial,*flags_final;
    float x,y,raper,raper1,gain=1.,phot_diam=6.,dis2,dis2_min,var_min=1.,flx0_fac=0.,spsf;
    float dx,dy,r,rmin,rmax,a,ai,la;
    int colnum_x,colnum_y,colnum_f,colnum_df,colnum_flags;
    int i,j,jj,ki,kj,xmin,xmax,ymin,ymax,*minx,*miny,*maxx,*maxy,*ixy,*norm,norm0,ki0,kj0,ki1,kj1;

    if (argc < 7) {
      printf("Usage: calc_phot datafile weightfile segfile catfile starfile psf_outfile [npsf] [crpix1] [crpix2] [refine_pos] [flx0_fac]\n");
      printf("    set npsf negative to skip calculation of psf flux (psf generation only)\n");
      printf("\n");
      return(0);
    }
    
    /* perhaps build the psf model */
    make_psf = access( argv[6], F_OK );
    /* perhaps calculate the psf model */
    if (argc>7) npsf = atoi(argv[7]);
    if (npsf<0) {
        make_calc=0;
        npsf = -npsf;
    }

    fits_open_file(&afptr, argv[1], READONLY, &status); /* open input image */
    fits_open_file(&bfptr, argv[2], READONLY, &status); /* open input weight image */
    fits_open_file(&cfptr, argv[3], READONLY, &status); /* open input image segmentation map */
    fits_open_file(&dfptr, argv[4], READWRITE, &status); /* open catalog table */
    fits_open_file(&efptr, argv[5], READONLY, &status); /* open psf star list */

    fits_get_img_dim(afptr, &anaxis, &status);  /* read dimensions */
    fits_get_img_size(afptr, 2, anaxes, &status);

    fits_get_img_size(efptr, 2, enaxes1, &status);
    npsf_stars = enaxes1[0];

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
      fits_read_key(afptr, TFLOAT, "PHOTDIAM", &phot_diam, NULL, &status);
      if (status!=0) {
          status=0;
          printf("Error reading phot_diam.\n");
          phot_diam=6.;
      }
      fits_read_key(afptr, TFLOAT, "CRPIX1", &x, NULL, &status);
      fits_read_key(afptr, TFLOAT, "CRPIX2", &y, NULL, &status);

      raper = phot_diam/2.;
      raper1 = 3*raper;
      irad = (int)ceil(2*raper);

      x0=0; y0=0;
      if (argc>8) x0=(int)(atof(argv[8])-x);
      if (argc>9) y0=(int)(atof(argv[9])-y);
      if (argc>10) refine_pos=atoi(argv[10]);
      if (argc>11) flx0_fac=atof(argv[11]);

      if (make_psf) {
          fits_read_key(efptr, TINT, "NX0", &nx0, NULL, &status);
          fits_read_key(efptr, TINT, "NY0", &ny0, NULL, &status);
          nbinx = npsf;
          nbiny = npsf;
      }
      else {
          fits_open_file(&outfptr, argv[6], READONLY, &status);
          fits_read_key(outfptr, TINT, "NX0", &nx0, NULL, &status);
          fits_read_key(outfptr, TINT, "NY0", &ny0, NULL, &status);
          fits_read_key(outfptr, TINT, "NBINX", &nbinx, NULL, &status);
          fits_read_key(outfptr, TINT, "NBINY", &nbiny, NULL, &status);
          fits_close_file(outfptr, &status);
      }
      if (status!=0) status=0;

      idx = (long *) malloc(npsf_stars * sizeof(long)); /* mem for psf star indices */
      if (idx==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }
      fits_read_pix(efptr, TLONG, firstpix, npsf_stars, NULL, idx, NULL, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */

      /* get table size */
      fits_movrel_hdu(dfptr, 1, NULL, &status);
      fits_get_num_rows(dfptr, &nrows, &status);

      /* find the columns numbers */
      fits_get_colnum(dfptr, CASEINSEN, "x_image", &colnum_x, &status);
      fits_get_colnum(dfptr, CASEINSEN, "y_image", &colnum_y, &status);
      fits_get_colnum(dfptr, CASEINSEN, "flux_aper", &colnum_f, &status);
      fits_get_colnum(dfptr, CASEINSEN, "fluxerr_aper", &colnum_df, &status);
      fits_get_colnum(dfptr, CASEINSEN, "flags", &colnum_flags, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */

      sz = 2*irad+1;
      nx = nbinx*sz;
      ny = nbiny*sz;

      mx = (float *) malloc(nrows * sizeof(float));
      my = (float *) malloc(nrows * sizeof(float));
      ixy = (int *) malloc(nrows * sizeof(int));
      flx = (float *) malloc(2*nrows * sizeof(float));
      dflx = (float *) malloc(2*nrows * sizeof(float));

      if (mx==NULL || my==NULL || flx==NULL || dflx==NULL || ixy==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      minx = (int *) malloc(nbinx*nbiny * sizeof(int));
      miny = (int *) malloc(nbinx*nbiny * sizeof(int));
      maxx = (int *) malloc(nbinx*nbiny * sizeof(int));
      maxy = (int *) malloc(nbinx*nbiny * sizeof(int));
      norm = (int *) malloc(nbinx*nbiny * sizeof(int));

      if (minx==NULL || miny==NULL || maxx==NULL || maxy==NULL || norm==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      nixy0 = (long *) malloc(nbinx*nbiny * sizeof(long));
      cnixy0 = (long *) malloc(nbinx*nbiny * sizeof(long));
      nixy_iter0 = (long *) malloc(nbinx*nbiny * sizeof(long));
      idx_ar0 = (long *) malloc(npsf_stars * sizeof(long)); /* mem for psf star indices */
      if (nixy0==NULL || cnixy0==NULL || nixy_iter0==NULL || idx_ar0==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      nixy = (long *) malloc(nbinx*nbiny * sizeof(long));
      cnixy = (long *) malloc(nbinx*nbiny * sizeof(long));
      nixy_iter = (long *) malloc(nbinx*nbiny * sizeof(long));
      idx_ar = (long *) malloc(nrows * sizeof(long)); /* mem for psf star indices */
      if (nixy==NULL || cnixy==NULL || nixy_iter==NULL || idx_ar==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      flags_initial = (int *) malloc(nrows * sizeof(int));
      flags_final = (int *) malloc(nrows * sizeof(int));
      psf = (float *) malloc(sz*sz * sizeof(float)); /* mem for image 1 */
      psf_stack = (float *) malloc(nx*ny * sizeof(float)); /* mem for image 1 */
     
      if (flags_initial==NULL || flags_final==NULL || psf==NULL || psf_stack==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      /* read in the table data */
      fits_read_col(dfptr, TFLOAT, colnum_x, 1, 1, nrows, NULL, mx, &anynulls, &status);
      fits_read_col(dfptr, TFLOAT, colnum_y, 1, 1, nrows, NULL, my, &anynulls, &status);
      fits_read_col(dfptr, TINT, colnum_flags, 1, 1, nrows, NULL, flags_initial, &anynulls, &status);
      fits_read_col(dfptr, TFLOAT, colnum_f, 1, 1, nrows*2, NULL, flx, &anynulls, &status);
      fits_read_col(dfptr, TFLOAT, colnum_df, 1, 1, nrows*2, NULL, dflx, &anynulls, &status);
      if (status) fits_report_error(stderr, status); /* print any error message */


      /* build ixy stamp locations and ordered list, idx_ar, of id's per stamp */
      for (kj=0;kj<nbiny;kj++) {
          for (ki=0;ki<nbinx;ki++) {
              jj = ki+nbinx*kj;
              maxx[jj]=maxy[jj]=0;
              minx[jj]=anaxes[0]; miny[jj]=anaxes[1];
              nixy0[jj]=nixy_iter0[jj]=nixy[jj]=nixy_iter[jj]=0;
          }
      }
      for (k=0;k<nrows;k++) {
          ixy[k] = (int)((mx[k]-1)*nbinx/nx0) + nbinx*(int)((my[k]-1)*nbiny/ny0);
          mx[k] -= x0; my[k] -= y0;
          jj = ixy[k];
          xmin = (int)(mx[k]-raper1-0.500001)-1;
          xmax = (int)(mx[k]+raper1+0.499999)+1;
          ymin = (int)(my[k]-raper1-0.500001)-1;
          ymax = (int)(my[k]+raper1+0.499999)+1;
          if (xmin<minx[jj]) minx[jj]=xmin;
          if (xmax>maxx[jj]) maxx[jj]=xmax;
          if (ymin<miny[jj]) miny[jj]=ymin;
          if (ymax>maxy[jj]) maxy[jj]=ymax;
          nixy[jj] += 1;
      }
      for (k0=0;k0<npsf_stars;k0++) {
          k=idx[k0];
          jj = ixy[k];
          nixy0[jj] += 1;
      }
      for (kj=0;kj<nbiny;kj++) {
          for (ki=0;ki<nbinx;ki++) {
              jj = ki+nbinx*kj;
              if (minx[jj]<0) minx[jj]=0;
              if (miny[jj]<0) miny[jj]=0;
              if (maxx[jj]>anaxes[0]-1) maxx[jj]=anaxes[0]-1;
              if (maxy[jj]>anaxes[1]-1) maxy[jj]=anaxes[1]-1;
              if (maxx[jj]-minx[jj]+1>npix_x) npix_x = maxx[jj]-minx[jj]+1;
              if (maxy[jj]-miny[jj]+1>npix_y) npix_y = maxy[jj]-miny[jj]+1;
          }
      }
      cnixy0[0]=cnixy[0]=0;
      for (jj=1;jj<nbinx*nbiny;jj++) {
          cnixy0[jj] = cnixy0[jj-1] + nixy0[jj-1];
          cnixy[jj] = cnixy[jj-1] + nixy[jj-1];
      }
      for (k=0;k<nrows;k++) {
          jj = ixy[k];
          mx[k] -= minx[jj];
          my[k] -= miny[jj];
          idx_ar[ nixy_iter[jj] + cnixy[jj] ] = k;
          nixy_iter[jj]+=1;
      }
      for (k0=0;k0<npsf_stars;k0++) {
          k = idx[k0];
          jj = ixy[k];
          idx_ar0[ nixy_iter0[jj] + cnixy0[jj] ] = k;
          nixy_iter0[jj]+=1;
      }

      npixels = npix_x*npix_y;  /* number of pixels per thumb */
      apix = (float *) malloc(npixels * sizeof(float));
      bpix = (float *) malloc(npixels * sizeof(float));
      cpix = (long *) malloc(npixels * sizeof(long)); /* segmentation map */

      if (apix == NULL || bpix==NULL || cpix==NULL) {
        printf("Memory allocation error\n");
        return(1);
      }

      if (refine_pos>0) {

          /* refine the source positions */
          fprintf(stderr," Refining source positions for %s...\n",argv[1]);

          for (kj=0;kj<nbiny;kj++) {
              for (ki=0;ki<nbinx;ki++) {
                  jj = ki+nbinx*kj;
                  if (maxx[jj]>minx[jj] && maxy[jj]>miny[jj]) {
                      firstpix[0] = minx[jj]+1; firstpix[1] = miny[jj]+1;
                      lastpix[0] = maxx[jj]+1; lastpix[1] = maxy[jj]+1;
                      fits_read_subset(afptr, TFLOAT, firstpix, lastpix, inc, NULL, apix, NULL, &status);
                  }
                  refine_mxy(maxx[jj]-minx[jj]+1,maxy[jj]-miny[jj]+1, nixy[jj], apix, idx_ar+cnixy[jj], mx, my, flx, dflx, phot_diam/3.);
              }
          }

      }

      if (make_psf) {

          /* get the array of psf's from the images */
          fprintf(stderr," Building PSFs for %s in %.2f pixel radius aperture...\n",argv[1],raper);

          norm0=0;
          for (kj=0;kj<nbiny;kj++) {
              for (ki=0;ki<nbinx;ki++) {
                  jj = ki+nbinx*kj;
                  if (maxx[jj]>minx[jj] && maxy[jj]>miny[jj]) {
                      firstpix[0] = minx[jj]+1; firstpix[1] = miny[jj]+1;
                      lastpix[0] = maxx[jj]+1; lastpix[1] = maxy[jj]+1;
                      fits_read_subset(afptr, TFLOAT, firstpix, lastpix, inc, NULL, apix, NULL, &status);
                      fits_read_subset(bfptr, TFLOAT, firstpix, lastpix, inc, NULL, bpix, NULL, &status);
                      fits_read_subset(cfptr, TLONG, firstpix, lastpix, inc, NULL, cpix, NULL, &status);
                      norm[jj] = build_median_psf(maxx[jj]-minx[jj]+1,maxy[jj]-miny[jj]+1, nixy0[jj], apix, bpix, cpix, psf, idx_ar0+cnixy0[jj], mx, my, flx, irad, gain*var_min);
                      for (j=0;j<sz;j++) {
                         j1 = ki*sz + sz*nbinx*(j + kj*sz);
                         for (i=0;i<sz;i++) psf_stack[i+j1] = psf[i+sz*j];
                      }
                  }
                  else norm[jj]=0;
                  if (norm[jj]<=0) norm0--;
              }
          }

          if (norm0<0) {
              /* fill in missing psf's with neighbors */
              for (kj=0;kj<nbiny;kj++) {
                  for (ki=0;ki<nbinx;ki++) {
                      if (norm[ki+nbinx*kj]<=0) {
                          dis2_min=1.e99;
                          for (kj1=0;kj1<nbiny;kj1++) {
                            for (ki1=0;ki1<nbinx;ki1++) {
                                dis2 = (kj-kj1)*(kj-kj1)+(ki-ki1)*(ki-ki1);
                                if (dis2<dis2_min && norm[ki1+nbinx*kj1]>0) {
                                  dis2_min=dis2;
                                  ki0=ki1; kj0=kj1;
                                }
                            }
                          }
                          if (dis2_min<1.e90) {
                            for (j=0;j<sz;j++) {
                               j1 = ki*sz + sz*nbinx*(j + kj*sz);
                               jj1 = ki0*sz + sz*nbinx*(j + kj0*sz);
                               for (i=0;i<sz;i++) psf_stack[i+j1] = psf_stack[i+jj1];
                            }
                          }
                      }
                  }
              }
          }
          for (k=0;k<=smooth_iter;k++) {
              for (kj=0;kj<nbiny-1;kj++) {
                  for (ki=0;ki<nbinx;ki++) {
                      for (j=0;j<sz;j++) {
                         j1 = ki*sz + sz*nbinx*(j + kj*sz);
                         jj1 = ki*sz + sz*nbinx*(j + (kj+1)*sz);
                         for (i=0;i<sz;i++) psf_stack[i+j1] = (psf_stack[i+j1]+psf_stack[i+jj1])/2.;
                      }
                      for (j=0;j<sz;j++) {
                         j1 = ki*sz + sz*nbinx*(j + (nbiny-1-kj)*sz);
                         jj1 = ki*sz + sz*nbinx*(j + (nbiny-2-kj)*sz);
                         for (i=0;i<sz;i++) psf_stack[i+j1] = (psf_stack[i+j1]+psf_stack[i+jj1])/2.;
                      }
                  }
              }
              for (kj=0;kj<nbiny;kj++) {
                  for (ki=0;ki<nbinx-1;ki++) {
                      for (j=0;j<sz;j++) {
                         j1 = ki*sz + sz*nbinx*(j + kj*sz);
                         jj1 = (ki+1)*sz + sz*nbinx*(j + kj*sz);
                         for (i=0;i<sz;i++) psf_stack[i+j1] = (psf_stack[i+j1]+psf_stack[i+jj1])/2.;
                      }
                      for (j=0;j<sz;j++) {
                         j1 = (nbinx-1-ki)*sz + sz*nbinx*(j + kj*sz);
                         jj1 = (nbinx-2-ki)*sz + sz*nbinx*(j + kj*sz);
                         for (i=0;i<sz;i++) psf_stack[i+j1] = (psf_stack[i+j1]+psf_stack[i+jj1])/2.;
                      }
                  }
              }
          }
      } else {
          fprintf(stderr," Reading existing psf file %s\n",argv[6]);
          fits_open_file(&outfptr, argv[6], READONLY, &status);
          firstpix[0]=firstpix[1]=1;
          fits_read_pix(outfptr, TFLOAT, firstpix, nx*ny, NULL, psf_stack, NULL, &status);
          fits_close_file(outfptr, &status);
      }
      if (make_psf) {
          /* save the psf model */
          fits_create_file(&outfptr, argv[6], &status);
          fits_copy_header(afptr, outfptr, &status);
          fits_update_key(outfptr, TINT, "NAXIS1", &nx, NULL, &status);
          fits_update_key(outfptr, TINT, "NAXIS2", &ny, NULL, &status);
          fits_update_key(outfptr, TINT, "NX0", &nx0, NULL, &status);
          fits_update_key(outfptr, TINT, "NY0", &ny0, NULL, &status);
          fits_update_key(outfptr, TINT, "NBINX", &nbinx, NULL, &status);
          fits_update_key(outfptr, TINT, "NBINY", &nbiny, NULL, &status);
          firstpix[0]=firstpix[1]=1;
          fits_write_pix(outfptr, TFLOAT, firstpix, nx*ny, psf_stack, &status);
          fits_close_file(outfptr, &status);
      } 

      if (make_calc) {

          fprintf(stderr," Fitting the PSFs to each source (%.2f)...\n",raper);

          for (kj=0;kj<nbiny;kj++) {
              for (ki=0;ki<nbinx;ki++) {
                  jj = ki+nbinx*kj;
                  if (maxx[jj]>minx[jj] && maxy[jj]>miny[jj]) {
                      spsf=0.;
                      for (j=0;j<sz;j++) {
                         j1 = ki*sz + sz*nbinx*(j + kj*sz);
                         dy = j-sz/2;
                         for (i=0;i<sz;i++) {
                             psf[i+sz*j] = psf_stack[i+j1];
                             dx = i-sz/2;
                             ai = MAX(fabs(dx),fabs(dy));
                             a = (ai>0) ? 0.5/ai : 0;
                             r = sqrt(dx*dx+dy*dy);
                             rmin = r*(1-a);
                             if (rmin<raper) {
                                 rmax = r*(1+a);
                                 la = (raper<rmax) ? (raper-rmin)/(rmax-rmin) : 1;
                                 spsf += la*psf[i+sz*j];
                             }
                         }
                      }
                      firstpix[0] = minx[jj]+1; firstpix[1] = miny[jj]+1;
                      lastpix[0] = maxx[jj]+1; lastpix[1] = maxy[jj]+1;
                      fits_read_subset(afptr, TFLOAT, firstpix, lastpix, inc, NULL, apix, NULL, &status);
                      fits_read_subset(bfptr, TFLOAT, firstpix, lastpix, inc, NULL, bpix, NULL, &status);
                      fits_read_subset(cfptr, TLONG, firstpix, lastpix, inc, NULL, cpix, NULL, &status);
                  }
                  calc_ap(maxx[jj]-minx[jj]+1,maxy[jj]-miny[jj]+1, nixy[jj], apix, bpix, cpix, idx_ar+cnixy[jj], mx, my, flx, dflx, raper, gain*var_min, flags_final);
                  calc_psf(maxx[jj]-minx[jj]+1,maxy[jj]-miny[jj]+1, nixy[jj], apix, bpix, cpix, psf, idx_ar+cnixy[jj], mx, my, flx, dflx, raper, gain*var_min, spsf);
                  calc_psf(maxx[jj]-minx[jj]+1,maxy[jj]-miny[jj]+1, nixy[jj], apix, bpix, cpix, psf, idx_ar+cnixy[jj], mx, my, flx, dflx, raper, gain*var_min, spsf);
                  calc_psf(maxx[jj]-minx[jj]+1,maxy[jj]-miny[jj]+1, nixy[jj], apix, bpix, cpix, psf, idx_ar+cnixy[jj], mx, my, flx, dflx, raper, gain*var_min, spsf);
              }
          }

          /* use flags_initial for deblending and flags_final for edges */
          for (k=0;k<nrows;k++) {
              if (flags_initial[k]>0 && flags_final[k]==0) flags_final[k] = 2;
          }

          /* now write out the new columns */
          fits_write_col(dfptr, TINT, colnum_flags, 1, 1, nrows, flags_final, &status);
          if (flx0_fac!=0) {
              /*printf (" Applying flx0_fac: %f\n",flx0_fac);*/
              for (k=0;k<nrows;k++) {
                  if (dflx[2*k]>0) flx[2*k] += flx0_fac*dflx[2*k];
              }
          }
          fits_write_col(dfptr, TFLOAT, colnum_f, 1, 1, nrows*2, flx, &status);
          fits_write_col(dfptr, TFLOAT, colnum_df, 1, 1, nrows*2, dflx, &status);

          if (refine_pos>0) {
              for (k=0;k<nrows;k++) {
                  jj = ixy[k];
                  mx[k] += minx[jj] + x0;
                  my[k] += miny[jj] + y0;
              }
              fits_write_col(dfptr, TFLOAT, colnum_x, 1, 1, nrows, mx, &status);
              fits_write_col(dfptr, TFLOAT, colnum_y, 1, 1, nrows, my, &status);
          }

          if (status) fits_report_error(stderr, status); /* print any error message */

      }

      free(apix);
      free(bpix);
      free(cpix);
      free(idx);
      free(mx);
      free(my);
      free(ixy);
      free(flx);
      free(dflx);

      free(nixy0);
      free(cnixy0);
      free(nixy_iter0);
      free(idx_ar0);
      free(nixy);
      free(cnixy);
      free(nixy_iter);
      free(idx_ar);

      free(minx);
      free(miny);
      free(maxx);
      free(maxy);
      free(norm);
      free(flags_initial);
      free(flags_final);
      free(psf);
      free(psf_stack);
    }
    fits_close_file(afptr, &status);
    fits_close_file(bfptr, &status);
    fits_close_file(cfptr, &status);
    fits_close_file(dfptr, &status);
    fits_close_file(efptr, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
