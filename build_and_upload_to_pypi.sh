#!/bin/bash

# Check if .env file exists
env_file=".env"

if [ ! -f "$env_file" ]; then
  echo "Error: .env file not found in the current directory."
  exit 1
fi

# Read and set environment variables from .env file
while IFS= read -r line || [[ -n "$line" ]]; do
  # Ignore lines starting with '#' (comments) or empty lines
  if [[ "$line" =~ ^\s*# || -z "$line" ]]; then
    continue
  fi

  # Split the line into key and value using the '=' delimiter
  key="${line%%=*}"
  value="${line#*=}"

  # Remove leading and trailing whitespace from the key and value
  key="${key// /}"
  value="${value// /}"

  # Set the environment variable
  export "$key=$value"
#   echo "Set environment variable: $key=$value"
done < "$env_file"

# # Upgrade pip
py -m pip install --upgrade pip

# Upgrade build
py -m pip install --upgrade build

# Build the project
py -m build

# Upgrade twine
py -m pip install --upgrade twine

# echo "USERNAME: $USERNAME"
# echo "API_TOKEN_TEST: $API_TOKEN_TEST"

# Upload the distribution files using twine

# Upload to test
twine upload --repository testpypi dist/*  -u "$USERNAME" -p "$API_TOKEN_TEST"

# Upload to PYPI prod
# twine upload dist/* -u "$USERNAME" -p "$API_TOKEN"

