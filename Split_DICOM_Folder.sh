#!/bin/bash

# Function to display script usage
usage() {
    echo "Usage: $0 -f /path/to/folder"
    exit 1
}

# Initialize variables
folder=""

# Parse command-line options
while getopts "f:" opt; do
    case $opt in
        f)
            folder="$OPTARG"
            ;;
        *)
            usage
            ;;
    esac
done

folder="${folder%/}"


# Check if the folder argument is provided
if [ -z "$folder" ]; then
    usage
fi

# Create the folder_2 directory next to the specified folder
folder_2="${folder}_2"

echo $folder_2
mkdir -p "$folder_2"

# Get a list of all files in the specified folder
files=("$folder"/*.MR.*.00??.????.??.??.*.IMA)

# Initialize variables to store the first two-digit number encountered
first_number=""

# Iterate over the files and extract the sequence number
for file in "${files[@]}"; do
    # Extract the sequence number using awk
	number=$(basename "$file" | awk -F'[._]' '{print $8}')
    # If it's the first number encountered, store it
    if [ -z "$first_number" ]; then
        first_number="$number"
    fi
 
    # Check if the current number is different from the first
    if [ "$number" != "$first_number" ]; then
        echo "Files in $folder have different sequence numbers."
    	sleep 1
        # Move files with a lower number to folder_2
        mv "$folder"/*.MR.*."$first_number".????.????.??.??.* "$folder_2"/
        echo "Files with sequence number $first_number moved to $folder_2."
        sleep 1
        # Exit the script
        exit 0
    fi
done

# If we reach this point, all files have the same two-digit number
echo "All files in $folder have the same sequence number ($first_number)."

