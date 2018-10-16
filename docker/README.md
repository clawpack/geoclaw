# Docker files

This folder contains docker files for creating images to test pull requests.

To create an image with name `test`, for example, with the docker file for 
[PR 8](https://github.com/barbagroup/geoclaw/pull/8):

```
$ docker build -t test -f Dockerfile-PR8 .
```

Then, to initialize a container with this image and get into an interactive Bash
shell, run:

```
$ docker run -it test
```

After getting into the shell, you can find the test case/example in 
`$CLAW/geoclaw/examples/landspill`. So you can run the example with
```
$ cd $CLAW/geoclaw/examples/landspill
$ make all
```
in the container.

### DockerHub image

Alternatively, you can also pull the image from DockerHub directly, instead of
building the image from the docker files.

For example, to use the image for [PR 8](https://github.com/barbagroup/geoclaw/pull/8):
from DockerHub, do:
```
$ docker pull barbagroup/landspill:PR8
```
and then create a container and get into the shell with
```
$ docker run -it barbagroup/landspill:PR8
```
