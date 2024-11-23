#!/bin/bash
# File: /data/proj2/home/students/m.borgmann/Master_thesis/systemscripts/project_summary.sh

# Create a timestamp for the filename
timestamp=$(date +"%Y%m%d_%H%M%S")
output_dir="$(dirname "$0")/outputs"
output_file="${output_dir}/project_summary_${timestamp}.md"

# Ensure the outputs directory exists
mkdir -p "$output_dir"

# Function to output to both file and console
output() {
    echo "$1" | tee -a "$output_file"
}

output "# Project Summary (Generated on $(date))"

output "## Project Structure:"
# List directories with depth of 5, excluding 'Python-3.9.9' and not expanding 'sra-toolkit'
tree /data/proj2/home/students/m.borgmann/Master_thesis -L 5 -I 'Python-3.9.9|sra-toolkit/*' |
    tee -a "$output_file"

output ""
output "## Python Version:"
python --version >> "$output_file" 2>&1

output "## Installed Packages:"
pip list >> "$output_file"

output "## Content of requirements.txt:"
cat "/data/proj2/home/students/m.borgmann/Master_thesis/requirements.txt" >> "$output_file"

output "## System Information:"
uname -a >> "$output_file"
echo "CPU: $(lscpu | grep 'Model name' | sed -r 's/Model name:\s{1,}//g')" >> "$output_file"
echo "Memory: $(free -h | awk '/^Mem:/ {print $2}')" >> "$output_file"

output "## Current Working Directory:"
pwd >> "$output_file"

output "## Active Virtual Environment:"
# Detect the active virtual environment (venv or conda)
if [[ -n "$VIRTUAL_ENV" ]]; then
    echo "$VIRTUAL_ENV" >> "$output_file"
elif [[ -n "$CONDA_DEFAULT_ENV" ]]; then
    conda info --envs | grep '*' >> "$output_file"
else
    echo "No active virtual environment" >> "$output_file"
fi

echo "Project summary has been saved to $output_file"
