ARG PYTHON_VERSION=3.11.2-slim-bullseye

FROM python:${PYTHON_VERSION}

RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates make gcc g++

WORKDIR /cleanser

COPY . .

RUN pip install cmdstanpy
RUN install_cmdstan
RUN pip install build
RUN python -m build .