#!/bin/bash
# Takes 

if [ $# -ne 5 ]; then
    echo "\
Usage: data-table0 [existing table] [I-band table] [I-band PCs]
                                    [V-band table] [V-band PCs]"
else
    table=$1
    tableI=$2
    PCI=$3
    tableV=$4
    PCV=$5

    # extract star names without file extensions from table
    sed 's/^\([^[:blank:]|.]*\).*$/\1/g' $tableI |
##    sed 's/^\([^ 	.]*\).[^ 	]*\([^\n]*\)/\1\2/g' $tableI |
##    sed '/.[a-zA-z]*/d' $tableI |
    paste - $PCI |



    # replace multiple spaces with single space
    sed 's/[[:blank:]]\{1,\}/ /g' |
    cat
fi
