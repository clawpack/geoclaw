This repository hosts the GeoClaw-GPU code described at: https://arxiv.org/abs/1901.06798.

Installation instructions for this GPU version:

1) Install GeoClaw by following the instructions here: http://www.clawpack.org/developers.html

2) For each of the submodules (amrclaw, geoclaw, clawutil, riemann, classic, pyclaw, visclaw), add a git remote that host corresponding GPU code:
    git remote add xinshengqin https://github.com/xinshengqin/amrclaw.git
    git remote add xinshengqin https://github.com/xinshengqin/clawutil.git
    git remote add xinshengqin https://github.com/xinshengqin/riemann.git
    git remote add xinshengqin https://github.com/xinshengqin/geoclaw.git
    git remote add xinshengqin https://github.com/xinshengqin/classic.git
    git remote add xinshengqin https://github.com/xinshengqin/pyclaw.git
    git remote add xinshengqin https://github.com/xinshengqin/visclaw.git

3) cd $CLAW and checkout following tags for each submodule:
    cd ./amrclaw  && git fetch xinshengqin && git checkout geo_gpu_paper && cd ../
    cd ./clawutil && git fetch xinshengqin && git checkout geo_gpu_paper && cd ../
    cd ./riemann  && git fetch xinshengqin && git checkout geo_gpu_paper && cd ../
    cd ./geoclaw  && git fetch xinshengqin && git checkout geo_gpu_paper && cd ../
    cd ./classic  && git fetch xinshengqin && git checkout geo_gpu_paper && cd ../
    cd ./pyclaw   && git fetch xinshengqin && git checkout geo_gpu_paper && cd ../
    cd ./visclaw  && git fetch xinshengqin && git checkout geo_gpu_paper && cd ../

4) Now one can run an example at geoclaw/examples/GPU/tsunami/chile2010/

Note that the code need PGI fortran compiler and CUDA (>=9.2) installed in the environment.

