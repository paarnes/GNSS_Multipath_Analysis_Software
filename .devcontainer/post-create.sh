#!/usr/bin/env bash
# Configure gh to use https. In Codespaces, gh is already authenticated via
# the injected GITHUB_TOKEN, so we only fall back to an interactive web login
# when no token is present (i.e. local Dev Containers).
set -e

gh config set git_protocol https >/dev/null || true

if [ -z "${GITHUB_TOKEN:-}" ] && [ -z "${GH_TOKEN:-}" ]; then
    if ! gh auth status >/dev/null 2>&1; then
        gh auth login --hostname github.com --git-protocol https --web || true
    fi
fi
