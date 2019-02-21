# Compile & run the code

Assume that LAPACK is installed correctly, just:

    cd Group
    make
    ./main.exe

For LAPACK you can either install/load it manually or use the pre-built Docker image.

# Use Docker image

See JiaweiZhuang/lapack_image/Dockerfile and https://cloud.docker.com/repository/docker/zhuangjw/lapack.

## On laptop

Install Docker and then:

    docker pull zhuangjw/lapack
    docker run --rm -it -v $(pwd):/work lapack
    cd /work

## On Harvard Odyssey

Use singularity. See more at https://www.rc.fas.harvard.edu/resources/documentation/software/singularity-on-odyssey/. Command are :

    srun -p test -n 1 -t 00-01:00 --pty --mem=4000 bash
    singularity pull docker://zhuangjw/lapack
    singularity shell lapack.simg
    alias ls="ls --color='auto'"

`make` inside singularity container will throw some warnings. Doesn't seem to affect anything.
