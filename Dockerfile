FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG JUST_VERSION=1.56.0
ARG POETRY_VERSION=2.2.1

RUN apt-get update \
    && apt-get install --yes --no-install-recommends \
        build-essential \
        ca-certificates \
        cmake \
        curl \
        git \
        gnupg \
        libmpfr-dev \
        parallel \
        python3 \
        python3-venv \
        software-properties-common \
        time \
        util-linux \
        zstd \
    && add-apt-repository --yes ppa:mhier/libboost-latest \
    && apt-get update \
    && apt-get install --yes --no-install-recommends libboost1.83-dev \
    && add-apt-repository --yes --remove ppa:mhier/libboost-latest \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 --silent --show-error --fail \
        https://just.systems/install.sh \
        | bash -s -- --tag "${JUST_VERSION}" --to /usr/local/bin

ENV CPM_SOURCE_CACHE=/opt/cpm-cache \
    MRNA_BUILD_JOBS=1 \
    MRNA_CLEAN_BUILDS=1 \
    MRNA_DIST=/opt/memerna \
    MRNA_TEST_COMPILERS=--compiler=default \
    PATH=/opt/memerna:/opt/poetry/bin:${PATH} \
    POETRY_HOME=/opt/poetry \
    POETRY_NO_INTERACTION=1 \
    POETRY_VIRTUALENVS_IN_PROJECT=1

RUN python3 -m venv "${POETRY_HOME}" \
    && "${POETRY_HOME}/bin/pip" install --no-cache-dir --quiet \
        "poetry==${POETRY_VERSION}" \
    && python3 --version \
    && poetry --version \
    && just --version

WORKDIR /src
COPY . .

RUN poetry install

RUN runtime_build="$(mktemp -d)" \
    && cmake -S . -B "${runtime_build}" -D CMAKE_BUILD_TYPE=Release \
    && cmake --build "${runtime_build}" --parallel "${MRNA_BUILD_JOBS}" \
        --target brute efn fold fuzz harness partition run_tests subopt \
    && "${runtime_build}/run_tests" \
    && mkdir -p "${MRNA_DIST}" \
    && install --mode=0755 \
        "${runtime_build}/brute" \
        "${runtime_build}/efn" \
        "${runtime_build}/fold" \
        "${runtime_build}/fuzz" \
        "${runtime_build}/harness" \
        "${runtime_build}/partition" \
        "${runtime_build}/subopt" \
        "${MRNA_DIST}/" \
    && cp --archive data "${MRNA_DIST}/data" \
    && rm -rf "${runtime_build}"

CMD ["just", "test"]
