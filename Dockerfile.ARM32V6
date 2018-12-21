#
# build our sysroot that we will cross compile against
#
FROM registry.gitlab.com/myriadrf/lime-suite:ARM32V6-master as sysroot

RUN apk update && \
	apk --force add \
		build-base \
		python3 \
		libusb-dev \
		fftw-dev \
		libstdc++ && \
	rm -rf /var/cache/apk/* && \
	pip3 install requests


#
# cross builder container for musl libc (aka alpine)
#
FROM registry.gitlab.com/pantacor/platform-tools/docker-musl-cross-arm as crossbuilder

WORKDIR /work
RUN mkdir /work/stage; apt-get update; apt-get install make cmake cmake-data -y; apt-get clean
COPY --from=sysroot / /sysroot-arm
COPY . src

RUN cd src && cmake --debug-output -DCMAKE_TOOLCHAIN_FILE=./cmake-cross/arm32v6 -DCMAKE_CXX_FLAGS="-I/sysroot-arm/usr/local/include" -DCMAKE_C_FLAGS="-I/sysroot-arm/usr/local/include"; make; make install

#
# produce our nice, tiny container with just the binary artifacts
#
FROM registry.gitlab.com/myriadrf/lime-suite:ARM32V6-master

COPY --from=crossbuilder /work/stage /usr/local

