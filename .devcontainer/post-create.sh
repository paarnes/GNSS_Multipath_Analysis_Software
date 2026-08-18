#!/usr/bin/env bash
set -e

# Editable install must happen here rather than in the Dockerfile, because the
# workspace is only bind-mounted after the image is built.
pip install --disable-pip-version-check --no-cache-dir -e .

# Configure gh to use https. In Codespaces, gh is already authenticated via
# the injected GITHUB_TOKEN, so we only fall back to an interactive web login
# when no token is present (i.e. local Dev Containers).
gh config set git_protocol https >/dev/null || true

if [ -z "${GITHUB_TOKEN:-}" ] && [ -z "${GH_TOKEN:-}" ]; then
    if ! gh auth status >/dev/null 2>&1; then
        gh auth login --hostname github.com --git-protocol https --web || true
    fi
fi
