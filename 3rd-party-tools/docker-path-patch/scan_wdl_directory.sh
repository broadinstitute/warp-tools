#!/bin/bash
# Scan a large set of WDL files with miniwdl in parallel. We only output a success line for WDLs with no errors. For WDLs
# with errors, we output an error line and the output of miniwdl.

ROOT_DIR=~/warp_patching/warp/

# Find all .wdl files in the $ROOT_DIR directory and its subdirectories
echo "Scanning files with miniwdl in parallel please wait..."
find "$ROOT_DIR" -name "*.wdl" -print0 | xargs -0 -n 1 -P "$(nproc)" -I {} sh -c '
    # Remove the leading directory $ROOT_DIR from the $file variable
    clean_file=${1#$2}

    # Run miniwdl on the file and redirect the output to /dev/null
    if miniwdl check --no-quant-check --suppress CommandShellCheck "$1" >/dev/null 2>&1
    then
        # If miniwdl returns a zero exit code, print "Success: {filename}"
        echo "Success: $clean_file"
    else
        # If miniwdl returns a non-zero exit code, print "Error: {filename}"
        echo "Error: $clean_file"
    fi
' sh {} "$ROOT_DIR" | tee >(grep "Success" > tmp__success.txt) >(grep "Error" > tmp__error.txt) >/dev/null

echo "Sorted Success Lines:"
sort -t ':' -k 1 tmp__success.txt
echo -e "\n\n-------------------------------------"
echo "Files with parsing errors:"
sort -t ':' -k 1 tmp__error.txt
echo -e "\n\n-------------------------------------"

# Run miniwdl on the errored files and print the output
while read -r file
do
    # Remove the "Error: " prefix from the filename
    clean_file=${file#"Error: "}

    # If miniwdl returns a non-zero exit code, print "Error: {filename}" and the output of miniwdl
    echo "Error: $clean_file"
    # miniwdl check --no-quant-check --debug --suppress CommandShellCheck "$ROOT_DIR$clean_file" 2>&1 | sed "s/^/   /"
    miniwdl check --no-quant-check --suppress CommandShellCheck "$ROOT_DIR$clean_file"
done < tmp__error.txt

rm tmp__success.txt tmp__error.txt
