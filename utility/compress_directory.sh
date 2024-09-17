#!/bin/bash

# Function to display help message
function show_help() {
    echo "Usage: $0 <directory> [options]"
    echo
    echo "Options:"
    echo "  -move                        Move the resulting .tar.gz archive to the current directory"
    echo "  -h, --help                   Show this help message and exit"
    echo
    echo "Description:"
    echo "  This script compresses the specified directory into a .tar.gz archive and deletes the original directory."
    echo "  By default, the archive will be kept in the same location as the original directory."
    echo "  If the -move option is specified, the archive will be moved to the current directory."
    echo
    echo "Example:"
    echo "  $0 /path/to/directory -move"
    echo
}

# Check if any argument is provided
if [[ $# -eq 0 ]]; then
    echo "Error: No directory provided."
    show_help
    exit 1
fi

# Parse command-line options
MOVE=false
while [[ "$1" != "" ]]; do
    case $1 in
        -h | --help )
            show_help
            exit 0
            ;;
        -move )
            MOVE=true
            ;;
        * )
            DIRECTORY=$1
            ;;
    esac
    shift
done

# Check if the provided directory exists
if [[ ! -d "$DIRECTORY" ]]; then
    echo "Error: Directory '$DIRECTORY' does not exist."
    exit 1
fi

# Get the directory path and basename
DIR_PATH=$(dirname "$DIRECTORY")
DIR_NAME=$(basename "$DIRECTORY")

# Compress the directory into a .tar.gz file in the same location as the original directory
ARCHIVE_NAME="$DIR_PATH/$DIR_NAME.tar.gz"
tar -czvf "$ARCHIVE_NAME" -C "$DIR_PATH" "$DIR_NAME"

# Check if compression was successful
if [[ $? -ne 0 ]]; then
    echo "Error: Failed to compress the directory."
    exit 1
fi

# Remove the original directory
rm -rf "$DIRECTORY"

# Move the archive if the -move option is specified
if [[ "$MOVE" == true ]]; then
    mv "$ARCHIVE_NAME" "$(pwd)"
    echo "Archive moved to $(pwd)/$DIR_NAME.tar.gz"
else
    echo "Archive created at $ARCHIVE_NAME"
fi

