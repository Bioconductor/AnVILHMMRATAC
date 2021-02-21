# Building the docker image

I built the [dockerhub image][] from the Dockerfile in this directory, using

```
docker build -t mtmorgan/hmmratac - < Dockerfile
```

The image is based on Alpine Linux. It installs tools necesary to
build bwa and samtools, as well as run Java software. It then adds the
source of each application, and compiles the source. The executables
are installed or linked under `/usr/local/bin` (so available on the
PATH) or (for the HMMRATAC.jar) in the login directory.

[dockerhub image]: https://hub.docker.com/repository/docker/mtmorgan/hmmratac
