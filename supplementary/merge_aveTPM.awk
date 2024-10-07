#!/bin/awk -f

BEGIN {
    FS = OFS = "\t";
    printf "annot";
}

# Process each input file
{
    # Extract sample name from the file name
    sampName = gensub(/\_aveTPM.*.txt$/, "", "1", FILENAME);

    # Store aveTPM value for each KO or EC
    aveTPM[sampName][$1] = $2;

    # Store descriptions
    descriptions[$1] = $3;
}

# Print header line with sample names
END {
    for (sampName in aveTPM) {
        printf "%s%s", OFS, sampName;
    }
    printf "%sDescription\n", OFS;

    # Print rows with aveTPM values
    for (i in descriptions) {
        printf "%s", i;
        for (sampName in aveTPM) {
            printf "%s%.5f", OFS, (i in aveTPM[sampName]) ? aveTPM[sampName][i] : 0;
        }
        printf "%s%s\n", OFS, descriptions[i];
    }
}
