#!/bin/sh

ra=244.002740071
dec=22.2690290469
radius=10.0

while [ $# -gt 0 ] ; do eval $1 ; shift ; done

[ "$ra" ] || ra=244.002740071
[ "$dec" ] || dec=22.2690290469

touch ps1_dr1_radec.txt

#url="http://vizier.u-strasbg.fr/cgi-bin/asu-tsv/?-source=II/349&-c.ra=${ra}&-c.dec=${dec}&-c.rm=${radius}&-out.max=10000"
url="https://vizier.u-strasbg.fr/viz-bin/asu-txt/?&-source=II/349&-c.ra=${ra}&-c.dec=${dec}&-c.rm=${radius}&-out.max=10000"
url="${url}&-out=RAJ2000&-out=DEJ2000&-out=gmag&-out=e_gmag&-out=rmag&-out=e_rmag&-out=imag&-out=e_imag&-out=zmag&-out=e_zmag&-out=ymag&-out=e_ymag&-sort=rmag"

wget $url -O ps1_dr1.txt

#grep -v '#' ps1_dr1.txt | awk '{if(NF==12) print}' | sed -n '4,$p' > ps1_dr1_radec.txt
[ -f ps1_dr1.txt ] && grep -v '#' ps1_dr1.txt | awk 'BEGIN{go=0}{if($1~/-------------/) go=1; if(NF>4 && go==1) print}' | sed -n '5,$p' > ps1_dr1_radec0.txt
if [ -f ps1_dr1_radec0.txt ]; then
    awk '{if ($0!~/[0-9]+\.[0-9]+\.[0-9]+/) print}' ps1_dr1_radec0.txt > ps1_dr1_radec0.tmp
    mv ps1_dr1_radec0.tmp ps1_dr1_radec0.txt
    awk '{if(NF==12) print}' ps1_dr1_radec0.txt > ps1_dr1_radec.txt
fi

if [ -s ps1_dr1_radec.txt ]; then

    awk '{print $1,$2,$3,$4}' ps1_dr1_radec.txt > ps1_dr1_radecg.txt
    awk '{print $1,$2,$5,$6}' ps1_dr1_radec.txt > ps1_dr1_radecr.txt
    awk '{print $1,$2,$7,$8}' ps1_dr1_radec.txt > ps1_dr1_radeci.txt
    awk '{print $1,$2,$9,$10}' ps1_dr1_radec.txt > ps1_dr1_radecz.txt
    awk '{print $1,$2,$11,$12}' ps1_dr1_radec.txt > ps1_dr1_radecy.txt
    # gri->BVRI from https://arxiv.org/pdf/1706.06147.pdf
    #
    C0=0.194; C1=0.561
    awk '{g=$3;dg=$4;r=$5;dr=$6; B = '$C0' + g + '$C1'*(g-r); dB = sqrt( ((1+'$C1')*dg)^2 + ('$C1'*dr)^2 ); print $1,$2,B,dB}' ps1_dr1_radec.txt > ps1_dr1_radecB.txt
    C0=-0.017; C1=-0.508
    awk '{g=$3;dg=$4;r=$5;dr=$6; V = '$C0' + g + '$C1'*(g-r); dV = sqrt( ((1+'$C1')*dg)^2 + ('$C1'*dr)^2 ); print $1,$2,V,dV}' ps1_dr1_radec.txt > ps1_dr1_radecV.txt
    C0=-0.166; C1=-0.275
    awk '{r=$5;dr=$6;i=$7;di=$8; R = '$C0' + r + '$C1'*(r-i); dR = sqrt( ((1+'$C1')*dr)^2 + ('$C1'*di)^2 ); print $1,$2,R,dR}' ps1_dr1_radec.txt > ps1_dr1_radecR.txt
    C0=-0.376; C1=-0.167
    awk '{g=$3;dg=$4;r=$5;dr=$6;i=$7;di=$8; I = '$C0' + i + '$C1'*(g-r); dI = sqrt( di^2 + ('$C1')^2*(dg^2+dr^2) ); print $1,$2,I,dI}' ps1_dr1_radec.txt > ps1_dr1_radecI.txt

    # now make some calibration bands for DDRAGO
    # B_D = B - 0.24*(B-V)
    # g_D = g + 0.05*(g-r)
    # r_D = r 
    # i_D = i + 0.04*(g-r)
    # gri_D = r + 0.07*(g-r)
    # z_D = z - 0.01*(r-i)
    # y_D = y 
    # zy_D = z - 0.09*(r-i)

    paste ps1_dr1_radecB.txt ps1_dr1_radecV.txt | awk '{print $1,$2,$3-0.24*($3-$7),$4}' > ps1_dr1_radecB_D.txt
    paste ps1_dr1_radecg.txt ps1_dr1_radecr.txt | awk '{print $1,$2,$3+0.05*($3-$7),$4}' > ps1_dr1_radecg_D.txt
    cp ps1_dr1_radecr.txt ps1_dr1_radecr_D.txt
    paste ps1_dr1_radeci.txt ps1_dr1_radecg.txt ps1_dr1_radecr.txt | awk '{print $1,$2,$3+0.04*($7-$11),$4}' > ps1_dr1_radeci_D.txt
    paste ps1_dr1_radecg.txt ps1_dr1_radecr.txt | awk '{print $1,$2,$7+0.07*($3-$7),$8}' > ps1_dr1_radecgri_D.txt
    paste ps1_dr1_radecz.txt ps1_dr1_radecr.txt ps1_dr1_radeci.txt | awk '{print $1,$2,$3-0.01*($7-$11),$4}' > ps1_dr1_radecz_D.txt
    cp ps1_dr1_radecy.txt ps1_dr1_radecy_D.txt
    paste ps1_dr1_radecz.txt ps1_dr1_radecr.txt ps1_dr1_radeci.txt | awk '{print $1,$2,$3-0.09*($7-$11),$4}' > ps1_dr1_radeczy_D.txt


else
    rm ps1_dr1_radec.txt 2>/dev/null
fi

rm ps1_dr1.txt 2>/dev/null
