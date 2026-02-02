#!/bin/bash

file=$1
shift

# get a dss file on the same pixel scale as input $file

dss_base_file=dss0.fits
setup_only=no

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

if [ -f "$dss_base_file" ]; then
    echo "Using existing dss_base_file $dss_base_file"
else

    ra0=`gethead CRVAL1 $file`
    dec0=`gethead CRVAL2 $file`
    cdec0=`echo $dec0 | awk '{printf("%.6f\n",cos($1/57.296))}'`
    x0=`gethead CRPIX1 $file`
    y0=`gethead CRPIX2 $file`
    dx=`gethead NAXIS1 $file`
    dy=`gethead NAXIS2 $file`
    cd11=`gethead CD1_1 $file`
    cd12=`gethead CD1_2 $file`
    cd21=`gethead CD2_1 $file`
    cd22=`gethead CD2_2 $file`

    gethead NAXIS1 NAXIS2 $file | awk '{\
        print '$ra0' + '$cd11'*(1-'$x0')+'$cd12'*(1-'$y0'),  '$dec0' + '$cd21'*(1-'$x0')+'$cd22'*(1-'$y0');\
        print '$ra0' + '$cd11'*(1-'$x0')+'$cd12'*($2-'$y0'),  '$dec0' + '$cd21'*(1-'$x0')+'$cd22'*($2-'$y0');\
        print '$ra0' + '$cd11'*($1-'$x0')+'$cd12'*(1-'$y0'),  '$dec0' + '$cd21'*($1-'$x0')+'$cd22'*(1-'$y0');\
        print '$ra0' + '$cd11'*($1-'$x0')+'$cd12'*($2-'$y0'),  '$dec0' + '$cd21'*($1-'$x0')+'$cd22'*($2-'$y0');\
    }' > radec$$.txt

    sort -n -k 1 radec$$.txt | awk '{if(NR==1) x0=$1; if(NR==4) {x1=$1; print 0.5*(x1+x0),(x1-x0)*60./'$cdec0'}}' > dra$$.txt
    sort -n -k 2 radec$$.txt | awk '{if(NR==1) x0=$2; if(NR==4) {x1=$2; print 0.5*(x1+x0),(x1-x0)*60}}' > ddec$$.txt

    ra1=`awk '{print $1}' dra$$.txt`
    dec1=`awk '{print $1}' ddec$$.txt`
    dra=`awk '{print $2}' dra$$.txt`
    ddec=`awk '{print $2}' ddec$$.txt`

    url="http://archive.stsci.edu/cgi-bin/dss_search?ra=${ra1}&dec=${dec1}&equinox=J2000&height=${ddec}&generation=2r&width=${dra}&format=FITS"
    echo "Grabbing dss from $url"
    wget $url -O $dss_base_file 1>/dev/null 2>&1

    rm dra$$.txt ddec$$.txt radec$$.txt

fi

if [ "$setup_only" = "yes" ]; then 
    echo "Setup performed, exiting."
    exit 3 
fi 

imhead -f $file > dss.head
swarp $dss_base_file -c ${SWARP_DIR}/coatli_redux.swarp -IMAGEOUT_NAME dss.fits -WEIGHT_TYPE NONE -WEIGHTOUT_NAME dss.wt.fits

base=${file%'.fits'}
wfile=${base}.wt.fits
if [ -f "$wfile" ]; then
    conpix -nr $wfile
    imarith dss.fits ${base}.wte.fits mul dss1.fits
    rm ${base}.wte.fits
    mv dss1.fits dss.fits
fi

rm dss.wt.fits dss.head
fits2jpg.py dss.fits 1024 linvert
