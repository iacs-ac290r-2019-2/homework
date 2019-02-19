If you want to install LAPACK/BLAS on Mac then just do it: https://pheiter.wordpress.com/2012/09/04/howto-installing-lapack-and-blas-on-mac-os/. Otherwise, Docker would be an easier solution.

(1) Install Docker

See https://docs.docker.com/docker-for-mac/install/

(2) Build image
    
    cd [directory containing Dockerfile]
    docker build -t lapack .

(3) Run image

    cd [directory containing your code]
    docker run --rm -it -v $(pwd):/work lapack
    cd /work

(4) Compile and run program

    gcc -o test_lapack.x test_lapack.c -llapacke
    ./test_lapack.x
