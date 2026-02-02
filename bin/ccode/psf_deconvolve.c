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

#define MIN(a,b) ((a) < (b) ? (a) : (b))

static double dsqrarg;
#define DSQR(a) ((dsqrarg=(a)) == 0.0 ? 0.0 : dsqrarg*dsqrarg)
#define SIGN(a,b) ((b) >= 0.0 ? fabs(a) : -fabs(a))

static double dmaxarg1,dmaxarg2;
#define DMAX(a,b) (dmaxarg1=(a),dmaxarg2=(b),(dmaxarg1) > (dmaxarg2) ? (dmaxarg1) : (dmaxarg2))

double dpythag(double a, double b) {
        double absa,absb;
        absa=fabs(a);
        absb=fabs(b);
        if (absa > absb) return absa*sqrt(1.0+DSQR(absb/absa));
        else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0+DSQR(absa/absb)));
}

void dsvdcmp(double *a, int n, double w[], double *v) {
	int flag,i,its,j,jj,k,l,nm;
	double anorm,c,f,g,h,s,scale,x,y,z,*rv1;

        rv1 = (double *) malloc(n * sizeof(double));
	g=scale=anorm=0.0;
	for (i=0;i<n;i++) {
		l=i+1;
		rv1[i]=scale*g;
		g=s=scale=0.0;
		if (i <= n) {
			for (k=i;k<n;k++) scale += fabs(a[k+n*i]);
			if (scale) {
				for (k=i;k<n;k++) {
					a[k+n*i] /= scale;
					s += a[k+n*i]*a[k+n*i];
				}
				f=a[i+n*i];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i+n*i]=f-g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=i;k<n;k++) s += a[k+n*i]*a[k+n*j];
					f=s/h;
					for (k=i;k<n;k++) a[k+n*j] += f*a[k+n*i];
				}
				for (k=i;k<n;k++) a[k+n*i] *= scale;
			}
		}
		w[i]=scale *g;
		g=s=scale=0.0;
		if (i <= n && i != n-1) {
			for (k=l;k<n;k++) scale += fabs(a[i+n*k]);
			if (scale) {
				for (k=l;k<n;k++) {
					a[i+n*k] /= scale;
					s += a[i+n*k]*a[i+n*k];
				}
				f=a[i+n*l];
				g = -SIGN(sqrt(s),f);
				h=f*g-s;
				a[i+n*l]=f-g;
				for (k=l;k<n;k++) rv1[k]=a[i+n*k]/h;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[j+n*k]*a[i+n*k];
					for (k=l;k<n;k++) a[j+n*k] += s*rv1[k];
				}
				for (k=l;k<n;k++) a[i+n*k] *= scale;
			}
		}
		anorm=DMAX(anorm,(fabs(w[i])+fabs(rv1[i])));
	}
	for (i=n-1;i>=0;i--) {
		if (i < n-1) {
			if (g) {
				for (j=l;j<n;j++) v[j+n*i]=(a[i+n*j]/a[i+n*l])/g;
				for (j=l;j<n;j++) {
					for (s=0.0,k=l;k<n;k++) s += a[i+n*k]*v[k+n*j];
					for (k=l;k<n;k++) v[k+n*j] += s*v[k+n*i];
				}
			}
			for (j=l;j<n;j++) v[i+n*j]=v[j+n*i]=0.0;
		}
		v[i+n*i]=1.0;
		g=rv1[i];
		l=i;
	}
	for (i=n-1;i>=0;i--) {
		l=i+1;
		g=w[i];
		for (j=l;j<n;j++) a[i+n*j]=0.0;
		if (g) {
			g=1.0/g;
			for (j=l;j<n;j++) {
				for (s=0.0,k=l;k<n;k++) s += a[k+n*i]*a[k+n*j];
				f=(s/a[i+n*i])*g;
				for (k=i;k<n;k++) a[k+n*j] += f*a[k+n*i];
			}
			for (j=i;j<n;j++) a[j+n*i] *= g;
		} else for (j=i;j<n;j++) a[j+n*i]=0.0;
		++a[i+n*i];
	}
	for (k=n-1;k>=0;k--) {
		for (its=1;its<=30;its++) {
			flag=1;
			for (l=k;l>=0;l--) {
				nm=l-1;
				if ((double)(fabs(rv1[l])+anorm) == anorm) {
					flag=0;
					break;
				}
				if ((double)(fabs(w[nm])+anorm) == anorm) break;
			}
			if (flag) {
				c=0.0;
				s=1.0;
				for (i=l;i<=k;i++) {
					f=s*rv1[i];
					rv1[i]=c*rv1[i];
					if ((double)(fabs(f)+anorm) == anorm) break;
					g=w[i];
					h=dpythag(f,g);
					w[i]=h;
					h=1.0/h;
					c=g*h;
					s = -f*h;
					for (j=0;j<n;j++) {
						y=a[j+n*nm];
						z=a[j+n*i];
						a[j+n*nm]=y*c+z*s;
						a[j+n*i]=z*c-y*s;
					}
				}
			}
			z=w[k];
			if (l == k) {
				if (z < 0.0) {
					w[k] = -z;
					for (j=0;j<n;j++) v[j+n*k] = -v[j+n*k];
				}
				break;
			}
			if (its == 30) fprintf(stderr,"no convergence in 30 dsvdcmp iterations\n");
			x=w[l];
			nm=k-1;
			y=w[nm];
			g=rv1[nm];
			h=rv1[k];
			f=((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
			g=dpythag(f,1.0);
			f=((x-z)*(x+z)+h*((y/(f+SIGN(g,f)))-h))/x;
			c=s=1.0;
			for (j=l;j<=nm;j++) {
				i=j+1;
				g=rv1[i];
				y=w[i];
				h=s*g;
				g=c*g;
				z=dpythag(f,h);
				rv1[j]=z;
				c=f/z;
				s=h/z;
				f=x*c+g*s;
				g = g*c-x*s;
				h=y*s;
				y *= c;
				for (jj=0;jj<n;jj++) {
					x=v[jj+n*j];
					z=v[jj+n*i];
					v[jj+n*j]=x*c+z*s;
					v[jj+n*i]=z*c-x*s;
				}
				z=dpythag(f,h);
				w[j]=z;
				if (z) {
					z=1.0/z;
					c=f*z;
					s=h*z;
				}
				f=c*g+s*y;
				x=c*y-s*g;
				for (jj=0;jj<n;jj++) {
					y=a[jj+n*j];
					z=a[jj+n*i];
					a[jj+n*j]=y*c+z*s;
					a[jj+n*i]=z*c-y*s;
				}
			}
			rv1[l]=0.0;
			rv1[k]=f;
			w[k]=x;
		}
	}
	free(rv1);
}

void dsvbksb(double *u, double w[], double *v, int n, double x[]) {
	int jj,j,i;
	double s,*tmp;

        tmp = (double *) malloc(n * sizeof(double));
	for (j=0;j<n;j++) {
		s=0.0;
		if (w[j]) {
			for (i=0;i<n;i++) s += u[i+n*j]*x[i];
			s /= w[j];
		}
		tmp[j]=s;
	}
	for (j=0;j<n;j++) {
		s=0.0;
		for (jj=0;jj<n;jj++) s += v[j+n*jj]*tmp[jj];
		x[j]=s;
	}
	free(tmp);
}

void basevec(int irad, double sg, int idx, int idy, double *fx, double *fy) {
    int i,dx,dy;
    double f0,s=1./(2*sg*sg);

    fx[irad] = (idx==0) ? 1 : 0;
    fy[irad] = (idy==0) ? 1 : 0;

    dx = (idx/2)*2 - idx;
    dy = (idy/2)*2 - idy;

    for (i=1;i<=irad;i++) {
       f0 = exp(-s*i*i);
       if (idx>0) {
           fx[irad+i] = f0*pow(i/sg,idx);
           fx[irad-i] = (dx==0) ? fx[irad+i] : -fx[irad+i];
       } else fx[irad+i]=fx[irad-i]=f0;
       if (idy>0) {
           fy[irad+i] = f0*pow(i/sg,idy);
           fy[irad-i] = (dy==0) ? fy[irad+i] : -fy[irad+i];
       } else fy[irad+i]=fy[irad-i]=f0;
    }
}

void convolve_vec(int irad, double *fx, double *fy, double *psf, double *cpsf0, double *cpsf) {
    int x,y,j,sz,j1,j2;
    long jj;
    double c0,fy0;

    sz=2*irad+1;

    for (y=0;y<sz;y++) {
        for (x=0;x<sz;x++) cpsf0[x+sz*y] = 0;
    }

    /* convolve in y */
    for (y=0;y<sz;y++) {
      for (j=-MIN(irad,2*irad-y);j<=MIN(y,irad);j++) {
          jj=sz*(y-j); fy0=fy[j+irad];
          for (x=0;x<sz;x++) cpsf0[x+sz*y] += psf[x+jj]*fy0;
      }
    }
    /* convolve in x */
    for (y=0;y<sz;y++) {
      for (x=0;x<sz;x++) {
          j1 = MIN(irad,2*irad-x); j2 = MIN(x,irad);
          c0=0;
          for (j=-j1;j<=j2;j++) c0 += cpsf0[x-j+sz*y]*fx[j+irad];
          cpsf[x+sz*y] = c0;
      }
    }
}

double dot_vec(int irad, double *psf, double *cpsf) {
    int i,j,sz=2*irad+1;
    double v0=0;
    for (j=0;j<sz;j++) {
        for (i=0;i<sz;i++) v0 += psf[i+j*sz]*cpsf[i+j*sz];
    }
    return v0;
}

int psf_deconv(int nn, int irad, double *fx, double *fy, double *psf_sci, double *psf_ref, double *mdl) {

   int i,j,k,l,sz,ssz;
   double *vec,*vec0,*vec1,*mat,*cpsf, *cpsf0;
   double *vm, *w, fy0, tmp;
   double lam=0.01,lam0=0.01,chi2=1.e90,chi0,chi2_last=1.e99,norm=1,norm_last=1,v1;

   sz=2*irad+1;
   ssz=sz*sz;

   vec = (double *) malloc(nn * sizeof(double));
   vec0 = (double *) malloc(nn * sizeof(double));
   vec1 = (double *) malloc(nn * sizeof(double));
   w = (double *) malloc(nn * sizeof(double));
   vm = (double *) malloc(nn*nn * sizeof(double));
   mat = (double *) malloc(nn*nn * sizeof(double));
   cpsf = (double *) malloc(nn*ssz * sizeof(double));
   cpsf0 = (double *) malloc(ssz * sizeof(double));

   chi0 = dot_vec(irad,psf_sci,psf_sci);

   for (k=0;k<nn;k++) {
       convolve_vec(irad,fx+sz*k,fy+sz*k,psf_ref,cpsf0,cpsf+ssz*k);
       vec[k] = dot_vec(irad,psf_sci,cpsf+ssz*k);
       for (l=0;l<=k;l++) mat[l+k*nn]=dot_vec(irad,cpsf+ssz*k,cpsf+ssz*l);
   }

   for (k=0;k<nn;k++) vec0[k] = sqrt(mat[k+k*nn]);
   for (k=0;k<nn;k++) {
       vec[k]/=vec0[k];
       mat[k+k*nn]=1+lam0;
       for (l=0;l<k;l++) mat[l+k*nn] /= vec0[k]*vec0[l];;
   }
   for (k=0;k<nn;k++) {
       for (l=0;l<k;l++) mat[k+l*nn] = mat[l+k*nn];
   }

   /* now solve mat*soln = vec */
   dsvdcmp(mat, nn, w, vm);

   for (k=0;k<nn;k++) {
     v1=0;
     for (l=0;l<nn;l++) v1 += mat[l+k*nn]*vec[l];
     vec1[k] = v1;
   }

   while (chi2*norm_last < chi2_last*norm && norm>0 && chi2>0 && lam>1.e-8) {
     chi2_last = chi2; norm_last=norm;
     chi2 = chi0; norm=1;
     for (k=0;k<nn;k++) {
       tmp=vec1[k]/(w[k]+lam-lam0);
       chi2-=tmp*tmp*(w[k]+2*lam-lam0);
       norm -= (w[k]-lam0)/(w[k]+lam-lam0)/nn;
     }
     if (norm>0 && chi2>0) {
       norm*=norm;
       lam/=2;
     }
   }

   for (k=0;k<nn;k++) w[k] += lam-lam0;
   dsvbksb(mat, w, vm, nn, vec);

   for (j=0;j<sz;j++) {
         for (i=0;i<sz;i++) mdl[i+sz*j]=0;
   }
   for (j=0;j<sz;j++) {
     for (k=0;k<nn;k++) {
         fy0 = fy[j+sz*k]*vec[k]/vec0[k];
         for (i=0;i<sz;i++) mdl[i+sz*j] += fx[i+sz*k]*fy0;
     }
   }

   free(w);
   free(vec);
   free(vec0);
   free(vec1);
   free(cpsf);
   free(cpsf0);
   free(vm);
   free(mat);

   return 0;
}


int main(int argc, char *argv[])
{
    fitsfile *psffile, *psffile0, *outfptr;
    int status = 0;  /* CFITSIO status value MUST be initialized to zero! */
    int i,j,k,l,k1,l1,irad, irad1, sz, sz1, nbinx, nbiny, nbinx1, nbiny1, nn=49,o,idx,idy;
    long pnaxes[2] = {1,1}, firstpix[2] = {1,1}, nx, ny, nx1, ny1;
    double *psf_stack, *delta_psf_stack_ref, *psf_stack_ref, *psf_ref, *psf_sci, *mdl, sg, *fx, *fy;
    float fwhm,fwhm0,sigma=1;
    const char fpre[] = "!delta_";
    char *ofile;

    if (argc < 3) {
      printf("Usage: psf_deconvolve psfgrid psfgrid_ref\n");
      printf("\n");
      return(0);
    }

    /* get the psf model */
    fits_open_file(&psffile, argv[1], READONLY, &status);
    fits_open_file(&psffile0, argv[2], READONLY, &status);
    fits_get_img_size(psffile0, 2, pnaxes, &status);
    nx = pnaxes[0]; ny = pnaxes[1];
    fits_get_img_size(psffile, 2, pnaxes, &status);
    nx1 = pnaxes[0]; ny1 = pnaxes[1];

    fits_read_key(psffile0, TINT, "NBINX", &nbinx, NULL, &status);
    fits_read_key(psffile0, TINT, "NBINY", &nbiny, NULL, &status);
    fits_read_key(psffile0, TFLOAT, "FWHM", &fwhm0, NULL, &status);
    fits_read_key(psffile, TINT, "NBINX", &nbinx1, NULL, &status);
    fits_read_key(psffile, TINT, "NBINY", &nbiny1, NULL, &status);
    fits_read_key(psffile, TFLOAT, "FWHM", &fwhm, NULL, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    status=0;

    sz = nx/(nbinx);
    irad = (sz-1)/2;
    sz1 = nx1/(nbinx1);
    irad1 = (sz1-1)/2;
    if (fwhm>fwhm0) sigma=sqrt(fwhm*fwhm-fwhm0*fwhm0)/2.3548;
    if (sigma<1) sigma=1;

    fprintf(stderr,"Convolving reference image psf to match science image psf with sigma=%.2f\n",sigma);

    psf_stack = (double *) malloc(nx1*ny1 * sizeof(double)); /* mem for psf stack */
    psf_stack_ref = (double *) malloc(nx*ny * sizeof(double)); /* mem for psf stack */
    delta_psf_stack_ref = (double *) malloc(nx*ny * sizeof(double)); /* mem for psf stack */
    mdl = (double *) malloc(sz*sz * sizeof(double));

    psf_sci = (double *) malloc(sz*sz * sizeof(double)); /* mem for psf */
    psf_ref = (double *) malloc(sz*sz * sizeof(double)); /* mem for psf */

    fx = (double *) malloc(nn*sz * sizeof(double));
    fy = (double *) malloc(nn*sz * sizeof(double));

    if (psf_stack==NULL || psf_stack_ref==NULL || delta_psf_stack_ref==NULL || psf_sci==NULL || psf_ref==NULL || mdl==NULL || fx==NULL || fy==NULL ) {
      printf("Memory allocation error\n");
      return(1);
    }

    fits_read_pix(psffile, TDOUBLE, firstpix, nx1*ny1, NULL, psf_stack, NULL, &status);
    fits_read_pix(psffile0, TDOUBLE, firstpix, nx*ny, NULL, psf_stack_ref, NULL, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */

    k=0;

    sg=0.5*sigma;
    for (o=0;o<=6;o++) {
        for (idx=0;idx<=o;idx++) {
            idy = o-idx;
            basevec(irad,sg,idx,idy,fx+k*sz,fy+k*sz);
            k++;
        }
    }
    sg=sigma;
    for (o=0;o<=4;o++) {
        for (idx=0;idx<=o;idx++) {
            idy = o-idx;
            basevec(irad,sg,idx,idy,fx+k*sz,fy+k*sz);
            k++;
        }
    }
    sg=2*sigma;
    for (o=0;o<=2;o++) {
        for (idx=0;idx<=o;idx++) {
            idy = o-idx;
            basevec(irad,sg,idx,idy,fx+k*sz,fy+k*sz);
            k++;
        }
    }

    /* loop over all psfs */
    for (l=0;l<nbiny;l++) {
      l1 = l*nbiny1/nbiny;
      for (k=0;k<nbinx;k++) {
        k1 = k*nbinx1/nbinx;

        for (j=-irad;j<=irad;j++) {
            for (i=-irad;i<=irad;i++) {
                if (abs(i)<=irad1 && abs(j)<=irad1) psf_sci[i+irad+sz*(j+irad)] = psf_stack[i+irad1+k1*sz1+nx1*(j+irad1+l1*sz1)];
                else psf_sci[i+irad+sz*(j+irad)] = 0;
                psf_ref[i+irad+sz*(j+irad)] = psf_stack_ref[i+irad+k*sz+nx*(j+irad+l*sz)];
            }
        }

        psf_deconv(nn,irad,fx,fy,psf_sci,psf_ref,mdl);
        for (j=0;j<sz;j++) {
            for (i=0;i<sz;i++) {
                psf_stack_ref[i+k*sz+nx*(j+l*sz)] = psf_ref[i+sz*j];
                delta_psf_stack_ref[i+k*sz+nx*(j+l*sz)] = mdl[i+sz*j];
            }
        }

      }
    }

    ofile = malloc(strlen(fpre)+strlen(argv[1])+1);
    strcpy(ofile,fpre); strcat(ofile,argv[1]);

    fits_create_file(&outfptr, ofile, &status);
    fits_copy_header(psffile, outfptr, &status);
    fits_update_key(outfptr, TINT, "NAXIS1", &nx, NULL, &status);
    fits_update_key(outfptr, TINT, "NAXIS2", &ny, NULL, &status);
    fits_update_key(outfptr, TINT, "NBINX", &nbinx, NULL, &status);
    fits_update_key(outfptr, TINT, "NBINY", &nbiny, NULL, &status);
    fits_write_pix(outfptr, TDOUBLE, firstpix, nx*ny, delta_psf_stack_ref, &status);
    fits_close_file(outfptr, &status);

    free(fx);
    free(fy);
    free(mdl);
    free(psf_stack);
    free(psf_stack_ref);
    free(delta_psf_stack_ref);
    free(psf_sci);
    free(psf_ref);
    free(ofile);

    fits_close_file(psffile, &status);
    fits_close_file(psffile0, &status);

    if (status) fits_report_error(stderr, status); /* print any error message */
    return(status);
}
