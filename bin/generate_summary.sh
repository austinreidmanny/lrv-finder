#!/bin/bash

# ---------------------------------------------------------------------------- #
# Table Summarizer
# ---------------------------------------------------------------------------- #
# This tool will take in a two-column tab delimited table of mapped reads per
# LRV species, and sum up the reads per species. It will then use those
# per-species counts to determine whether or not each LRV is present in the run
# ---------------------------------------------------------------------------- #

# Make the required table is provided
if [[ "$#" < 1 ]]; then
    echo ""
    echo "WARNING: No table provided for analysis."
    echo "Usage: $0 initial_table.txt"
    echo "Retry with a table file provided..."
    echo ""
    exit 1
fi

# Read in the input table
input=$1

# Set up an intermediate and an output table
temp_table="temp_table.txt"
output_table="final_results.txt"

# Establish the "presence threshold" -- how many reads does it take to mean that the virus is present?
presence_threshold=10

# Read through read counts for each virus strain and sum up read counts per species
cat $input |
awk -v ONE=0 -v TWO=0 -v THREE=0 -v FOUR=0 -v FIVE=0 -v TV=0 -v TT=0 '{

    if ($1 ~ "LRV1")
        ONE += $2;
    else if ($1 ~ "LRV2")
        TWO += $2;
    }

END {
    print "LRV1" "\t" ONE "\n"       \
          "LRV2" "\t" TWO "\n"       \
    }' > "$temp_table"

# Make a header for the output table
echo -e "\n-------------------------------------------------------" > "$output_table"
echo "Raw read counts for Leishmaniavirus spp. in this run" >> "$output_table"
echo " (species with >= ${presence_threshold} reads are considered 'present')" >> "$output_table"
echo "-------------------------------------------------------" >> "$output_table"

# Copy intermediate table to final table; add a blank line at the end
cat "$temp_table" >> "$output_table"
echo -e "-------------------------------------------------------\n" >> "$output_table"
echo -e "-------------------------------------------------------" >> "$output_table"
echo "Summary of Leishmaniavirus spp. present in this run" >> "$output_table"
echo "-------------------------------------------------------" >> "$output_table"

# Determine whether that final table contains enough reads for each species
while read -r species counts; do

    # Skip line if it doesn't have a counts value
    if [[ -z "$counts" ]]; then
        echo ""

    elif [[ $counts -gt $presence_threshold ]]; then
        echo "$species present"

    else
        echo "$species not present"
    fi

done >> "$output_table" \
< "$temp_table"

# Add one final formatting line at the end
echo "-------------------------------------------------------" >> "$output_table"

# Remove intermediate table
rm "$temp_table"
