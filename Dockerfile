ARG PYTHON_VERSION=3.11.2-slim-bullseye

FROM gcc:13.2.0-bookworm as stan

ENV STAN_VER=2.33.1
ENV CMDSTAN=cmdstan-${STAN_VER}

WORKDIR /
RUN apt-get update && apt-get install --no-install-recommends -qq wget ca-certificates

RUN wget https://github.com/stan-dev/cmdstan/releases/download/v${STAN_VER}/${CMDSTAN}.tar.gz
RUN tar -zxpf ${CMDSTAN}.tar.gz
RUN mv ${CMDSTAN} cmdstan

WORKDIR /cmdstan
RUN make build

COPY ./guide-mixture.stan .
RUN make guide-mixture

FROM python:${PYTHON_VERSION}

COPY --from=stan /cmdstan/guide-mixture /usr/local/bin

WORKDIR /cleanser

COPY . .

RUN pip install build
RUN python -m build
