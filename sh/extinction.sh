#!/bin/bash

if [ $# -ne 3 ]
then
    echo \
"Usage: extinction [ident.dat] [reddening map] [star type]

ident.dat file taken from OGLE
reddening map must be in FITS format
star type must be the type identifier given in the 4th column of ident.dat"
    exit 1
fi

ident=$1
map=$2
type=$3

stars=$(mktemp)
radec=$(mktemp)
xy=$(mktemp)
ext=$(mktemp)
awk -v type=$type '$4 == type' $ident > $stars
awk '{print $5, $6}' $stars > $radec
sky2xy $map @$radec | awk '{print $5,$6}' > $xy
getpix -v $map @$xy | awk '{print $3}' > $ext

awk '{print $1}' $stars |
paste - $ext |
sed 's/[[:blank:]]\{1,\}/ /g'

rm $stars $radec $xy $ext
