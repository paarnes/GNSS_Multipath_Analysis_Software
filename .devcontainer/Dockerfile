# Debian 13 "trixie" based Python image.
# Python 3.14 is not yet published as an official devcontainer variant,
# so we pin to the latest available trixie variant (3.13).
ARG VARIANT="3.13-trixie"
FROM mcr.microsoft.com/devcontainers/python:${VARIANT}

# Install GitHub CLI
RUN apt-get update && export DEBIAN_FRONTEND=noninteractive \
    && apt-get -y install --no-install-recommends \
        curl \
        gnupg \
    && curl -fsSL https://cli.github.com/packages/githubcli-archive-keyring.gpg | \
        gpg --dearmor -o /usr/share/keyrings/githubcli-archive-keyring.gpg \
    && echo "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/githubcli-archive-keyring.gpg] https://cli.github.com/packages stable main" \
        | tee /etc/apt/sources.list.d/github-cli.list > /dev/null \
    && apt-get update \
    && apt-get install -y gh \
    && rm -rf /var/lib/apt/lists/*

# Python tooling: upgrade pip, install project requirements + Jupyter
COPY .devcontainer/requirements.txt /tmp/pip-tmp/
RUN python -m pip install --upgrade pip \
    && pip3 --disable-pip-version-check --no-cache-dir install -r /tmp/pip-tmp/requirements.txt \
    && pip3 --disable-pip-version-check --no-cache-dir install notebook jupyterlab ipykernel \
    && rm -rf /tmp/pip-tmp
