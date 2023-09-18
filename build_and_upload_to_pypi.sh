#!/bin/bash

# Source the .env file to load environment variables
if [ -f .env ]; then
  source .env
else
  echo "Error: .env file not found."
  exit 1
fi

# Ensure that USERNAME and API_TOKEN are set
if [ -z "$USERNAME" ] || [ -z "$API_TOKEN" ]; then
  echo "Error: USERNAME or API_TOKEN is not set in .env."
  exit 1
fi

# Upgrade pip
py -m pip install --upgrade pip

# Upgrade build
py -m pip install --upgrade build

# Build the project
py -m build

# Upgrade twine
py -m pip install --upgrade twine

# Upload the distribution files using twine
TWINE_USERNAME="$USERNAME" TWINE_PASSWORD="$API_TOKEN_TEST" twine --repository testpypi dist/*
# TWINE_USERNAME="$USERNAME" TWINE_PASSWORD="$API_TOKEN" twine upload dist/*
