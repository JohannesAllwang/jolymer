#! /bin/bash

shopt -s -o nounset

#
# Smoothing and rebinning a SAXS curve to use it in SAXS-driven MD.
#
# Explantion: In SAXS-driven MD, using more than 1,5 or 2 q-points per Shannel channel has no noticible effect
#             except that the simulations would run drastically slower.
#
# Input: experimental, noisy and oversampled SAXS curve with typically many data points.
# Ouput: Smoothed SAXS curve, with typically only 2 q-points per Shannon channel. The new error bars
#        reflects the uncertainty of I(q) AFTER averaging over neihboring q-points within each Shannon bin.
#
# Written by Milos Ivanovic and Jochen Hub, Saarland University, 2021
#
#
# Licence: GNU General Public License, version 3 or later
#
# This program is free software: you can redistribute it and/or modify it under the terms of the
# GNU General Public License as published by the Free Software Foundation, either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without
# even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details, https://www.gnu.org/licenses/ 
#

# pick column(s) of data file or line
column ()
{
    awk -v cols="$1" '{split(cols,c,","); for (i=1; i<length(c); i++) {printf "%s ", $c[i]}; print $c[length(c)] }'
}

# command line calculator
= ()
{
    local in="$(echo "$@" | sed -e 's/\[/(/g' -e 's/\]/)/g')";
    awk 'BEGIN {print '"$in"'}' < /dev/null
}

# ATSAS in PATH ?
for tool in autorg datgnom datregrid shanum; do
    if ! which $tool >& /dev/null; then
        echo "Error, $tool not found in path. Install the ATSAS software and add it to your PATH variable."
        exit 1
    fi
done


f=""                       # input SAXS curve
o=saxs_rebinned.xvg        # output rebinned SAXS curve
m=2                        # Number of q-points written per Shannon bin. For a target curve for SAXS-driven MD, m=2 is enough.
osm=smoothed_allpoints.xvg # smoothed curve, but still all points (optional)
bVerbose=no
bSkipGnom=no               # only rebinning, without smothing with GNOM first. This way, the rebinned SAXS curve gets a bit "bumpy".
                           # Does not work in case that SHANUM fails, as ocasionally happening

while [ $# -gt 0 ]; do
    case "$1" in
        -f)          shift;            f="$1"  ;;
        -osm)        shift;          osm="$1"  ;;
        -o)          shift;            o="$1"  ;;
        -m)          shift;            m=$1    ;;
        -skip-gnom)            bSkipGnom=yes   ;;
        -v)                     bVerbose=yes   ;;
        
        *)           echo -e "\nError, unknown argument: $1"; exit 1  ;;
    esac
    shift
done

if [ "$f" = "" ]; then
    echo -e "Usage:\nsmooth-saxs-curve.sh -f saxs.dat -o [saxs_rebinned.dat|saxs_rebinned.xvg] [-osm smoothed_allpoints.dat] [-m POINTS_PER_SHANNON_BIN]"; exit 1
fi
if [ ! -e "$f" ]; then
    echo -e "\nError, input SAXS curve $f not found\n"; exit 1
fi

rm -f temp_in out_shanum gnom.out autorg.out $o datgnom.out
echo "NOTE: Assuming that q points in $f are equidistant."

egrep -v '#|@|&' $f > temp_in

if [ $bSkipGnom = no ]; then
    #
    # Get Rg with AUTORG, needed for DATGNOM
    #
    autorg temp_in >& autorg.out || exit 1
    
    if [ $bVerbose = yes ]; then
        echo -e "AUTORG output:\n****************"
        cat autorg.out
        echo "****************"
    fi
    
    Rg=$(grep ^Rg autorg.out | column 3)
    
    #
    # Now smoothing with DATGNOM
    #
    echo "Getting smoothed curve with DATGNOM ..."
    datgnom $f -rg $Rg -o datgnom.out --last 1000000
    Dmax=$(grep "Maximum characteristic size" datgnom.out | column 4)

    
    # echo -e "\n${f}\n\n\n\n${Dmax}\n\n\n\ngnom.out" | gnom >& gnom.err
    sed -n '/J EXP/,/Real Space/ p' < datgnom.out | egrep -v 'EXP|Dist|Data' | grep [0-9] | awk 'NF==5' | awk '{print $1, $4,$3}' > $osm

    extension=${osm##*.}
    if [[ $extension = xvg || $extension = agr ]]; then
        {
            echo "@type xydy"
            cat $osm
        } > tmp123 && mv tmp123 $osm
    fi
    
    cp $osm temp_in
fi

echo "Picking up number of Shannon channels with SHANUM ..."
shanum temp_in &> out_shanum

if grep -q Nopt out_shanum; then
    Nq=`grep Nopt out_shanum | column 2`
else
    qmin=$(grep [0-9] temp_in | head -n1 | column 1)
    qmax=$(grep [0-9] temp_in | tail -n1 | column 1)
    Nq=$(= "($qmax-$qmin)*$Dmax/3.14159")
    printf "\nNOTE: SHANUM could not detect the number of Shannon channels.\n"
    printf   "      Computing number of Shannon channels instead via (q[max]-q[min])*D/pi = %g\n\n" $Nq
fi


#
# Join $ neighboring I(q) points, and compute sigma after noining $njoin points
#
nq_total=`wc -l < temp_in`
q_Shan=$(= "$nq_total/$Nq")
njoin=$(= "$q_Shan/$m" | awk '{printf "%.0f\n", $1}')

echo "Rebinning with DATREGRID ..."
datregrid --join=${njoin} temp_in --output=datregrid_out
grep -v arent datregrid_out | egrep -v 'creator|datregrid-operation|Sample'  > $o

#
# In our SAXS-driven MD, sigma must represent the error after averaging within one Shannon bin. In case we
# write more than one data point per Shannon bin ($m>1), the error written by datregrid
# is too large by sqrt($m), so correct the error:
#
echo "Error in $o is the error after averaging within one Shannon bin ..."
awk -v m=$m '{printf "%10.6e  %10.6e  %10.6e\n", $1, $2, $3/sqrt(m)}' < $o > tmp123 && mv tmp123 $o

#
# Add xydy key to xmgrace files
#
extension=${o##*.}
if [[ $extension = xvg || $extension = agr ]]; then
    {
        echo "@type xydy"
        cat $o
    } > tmp123 && mv tmp123 $o
fi

printf "\n\n"
printf "Rg   by AUTORG                                             = %g\n" $Rg
printf "Dmax by DATGNOM                                            = %g\n" $Dmax
printf "Number of q points read                                    = %d\n" $nq_total
printf "Number of independent data points according to SHANUM      = %g\n" $Nq
printf "Read    this number of q-points per Shannon channel        = %g\n" $q_Shan
printf "Writing this number of q-points per Shannon channel (-m)   = %g\n" $m
printf "Averaging over this number of q points                     = %d\n" $njoin

printf "\nOutput:\n-------\n"
printf "Smoothed curve, still with all q-points       :  $osm\n"
printf "Smoothed curve, rebinned (for SAXS-driven MD) :  $o\n\n"

# Clean up
if [ $bVerbose = yes ]; then
    echo "Verbose, keeping all temporary files."
else
    rm -f  datregrid_out autorg.out temp_in out_shanum gnom.out temp_in
fi
