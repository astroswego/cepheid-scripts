#!/bin/bash
# Takes an existing data table, where each row of data has in its first column
# a unique ID for the row, and a series of files of names and data to insert
# into the specified columns

if (($# && !(($#+1) % 3)))
then
    # file for original data table
    origfile=$1
    # number of columns used in the original data table for the field
    # which is being updated. 0 if the field is new
    existingcols=$2

    # make temporary files for storing intermediate results
    infile=$(mktemp)
    outfile=$(mktemp)
    
    # remove heading from original file
    tail -n +2 $origfile > $outfile
    # count number of columns in original file
    origcols=$(head -n 2 $origfile | tail -n 1 | wc -w)
    # highest column that has been altered
    highestcol=0
    # amount of extra columns added so far
    colshift=0

    # iterate over each set of 3 arguments, where the first is a file
    # of star names, the second is a file of data associated with those
    # star names, and the third is the first column to insert the data at
    for (( i = 3; i <= $#; i += 3 ))
    do
        eval "namefile=\$$i
              datafile=\$$((i+1))
              tablecol=\$$((i+2))"
	if [ $tablecol -le $highestcol ]
	then
	    echo "Arguments must be given in ascending order."
	    # delete tmp files
	    rm $infile $outfile
	    exit 1
	else
	    # record the highest column used in the original table
	    highestcol=$tablecol
	fi
	# count columns in data file
	datacols=$(head -n 1 $datafile | wc -w)
	# copy output file to input file so that it can be used for output
	cp $outfile $infile
	# extract names from namefile, removing any file extensions
	sed "s/[.[:blank:]].*//" $namefile |
	# paste names with their data
	paste - $datafile |
	join $infile - > $outfile
	# output file becomes new input file
	cp $outfile $infile
	# now use awk to overwrite the appropriate columns, and save the
	# result to the output file
	awk -v existingcols=$existingcols \
            -v tablewidth=$((origcols+colshift)) \
            -v startcol=$((tablecol+colshift)) \
            -v newcols=$datacols '{
        for (col = 1;                     col <  startcol;           col++)
            printf $col " ";
        for (col = tablewidth+1;          col <= tablewidth+newcols; col++)
            printf $col " ";
        for (col = startcol+existingcols; col <= tablewidth;         col++)
            printf $col " ";
        printf "\n"
        }' $infile > $outfile
	# account for the change in the number of columns in the table
	((colshift += datacols-existingcols))
    done
    # replace all whitespace with tabs
    # repeated whitespace is replaced with a single tab
    # the result of this is what is output to stdout
    sed 's/[[:blank:]]\{1,\}/	/g' $outfile
    # delete tmp files
    rm $infile $outfile
else
    echo "\
Usage: data-table0 [existing table] [existing PC #]
                   ([name table] [data table] [column #])*

  column #'s must occur in ascending order."
fi
