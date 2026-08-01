FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive
ARG JUST_VERSION=1.56.0
ARG UV_VERSION=0.11.29

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
        software-properties-common \
        time \
        tmux \
        util-linux \
        zstd \
    && add-apt-repository --yes ppa:mhier/libboost-latest \
    && apt-get update \
    && apt-get install --yes --no-install-recommends libboost1.83-dev \
    && add-apt-repository --yes --remove ppa:mhier/libboost-latest \
    && rm -rf /var/lib/apt/lists/*

RUN curl --proto '=https' --tlsv1.2 --location --silent --show-error --fail \
        "https://astral.sh/uv/${UV_VERSION}/install.sh" \
        | env UV_UNMANAGED_INSTALL=/usr/local/bin sh
RUN curl --proto '=https' --tlsv1.2 --silent --show-error --fail \
        https://just.systems/install.sh \
        | bash -s -- --tag "${JUST_VERSION}" --to /usr/local/bin
RUN uv --version && just --version

ENV CPM_SOURCE_CACHE=/opt/cpm-cache \
    MRNA_DIST=/opt/memerna \
    MRNA_BUILD_JOBS=2 \
    MRNA_CLEAN_BUILDS=1 \
    MRNA_TEST_COMPILERS=--compiler=default \
    PATH=/opt/memerna:${PATH} \
    UV_LINK_MODE=copy \
    UV_PROJECT_ENVIRONMENT=/opt/venv \
    UV_PYTHON_INSTALL_DIR=/opt/uv/python

WORKDIR /src
COPY . .

RUN uv sync --locked

RUN runtime_build="$(mktemp -d)" \
    && cmake -S . -B "${runtime_build}" -D CMAKE_BUILD_TYPE=Release \
    && cmake --build "${runtime_build}" --parallel "${MRNA_BUILD_JOBS}" \
        --target efn fold fuzz harness partition run_tests subopt \
    && "${runtime_build}/run_tests" \
    && mkdir -p "${MRNA_DIST}" \
    && install --mode=0755 \
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
