FROM gcc:14.2

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    ninja-build \
    pkg-config \
    zlib1g-dev \
    libbz2-dev \
    liblzma-dev \
    libcurl4-openssl-dev \
    libdeflate-dev \
    wget \
    && rm -rf /var/lib/apt/lists/*

# Install htslib
RUN wget -q https://github.com/samtools/htslib/releases/download/1.21/htslib-1.21.tar.bz2 \
    && tar -xjf htslib-1.21.tar.bz2 \
    && cd htslib-1.21 && ./configure && make -j$(nproc) && make install \
    && cd / && rm -rf htslib-1.21 htslib-1.21.tar.bz2 \
    && ldconfig

WORKDIR /src

# Copy source
COPY CMakeLists.txt .
COPY include/ include/
COPY src/ src/

# Build
RUN cmake -B build -S . -G Ninja -DCMAKE_BUILD_TYPE=Release -DBUILD_TESTING=OFF \
    && cmake --build build -j$(nproc)

ENTRYPOINT ["/src/build/atroplex"]