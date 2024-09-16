#!/bin/bash

# Function to display help message
function show_help() {
    echo "Usage: $0 [options]"
    echo
    echo "Options:"
    echo "  -h, -help   Show this help message and exit"
    echo
    echo "Description:"
    echo "  This script searches for symbolic links in the current directory and all subdirectories."
    echo "  It displays symbolic links along with details such as user, group, permissions, and time of creation."
    echo
}

# Check if -h or -help is passed
if [[ "$1" == "-h" || "$1" == "-help" ]]; then
    show_help
    exit 0
fi

# Search for symbolic links and show detailed information
echo "Searching for symbolic links:"

# Find all symbolic links (-type l) and display detailed information with ls -l (permissions, user, group, time, target)
find . -type l -exec ls -l {} +

# Print message if no symbolic links are found
if [ $? -ne 0 ]; then
    echo "No symbolic links found."
fi

