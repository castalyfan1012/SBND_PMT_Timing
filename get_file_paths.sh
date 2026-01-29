#!/bin/bash

# Define the query range
RUN_RANGE="19403-19467"

# Output file
OUTPUT_FILE="file_list.txt"

# Clear the output file first
> "$OUTPUT_FILE"

echo "Listing files for run range $RUN_RANGE..."
FILES=$(samweb -e sbnd list-files "file_type data and data_tier raw and run_number $RUN_RANGE")

echo "Found $(echo "$FILES" | wc -l) files."
echo "Locating each file path..."

# Loop over each file and locate its path
while read -r FILE; do
    if [[ -n "$FILE" ]]; then
        LOCATION=$(samweb locate-file "$FILE" 2>/dev/null | grep '^enstore:' | head -n1)
        # Clean up the output (remove enstore: prefix and parentheses)
        CLEAN_PATH=$(echo "$LOCATION" | sed -E 's/^enstore://; s/\(.*\)//; s/ *$//')
        echo "$CLEAN_PATH/$FILE" >> "$OUTPUT_FILE"
        echo "✓ $FILE"
    fi
done <<< "$FILES"

echo
echo "Done. Paths saved to: $OUTPUT_FILE"
