#!/usr/bin/env bash

# Configure CI test environment variables for GeoClaw workflows.
# Usage:
#   source .github/scripts/set_test_env.sh regression <compiler> <build>

set -euo pipefail

test_group="${1:-}"
compiler="${2:-}"
build="${3:-}"

if [[ -z "${test_group}" || -z "${compiler}" || -z "${build}" ]]; then
  echo "Usage: source .github/scripts/set_test_env.sh <test_group> <compiler> <build>"
  return 1 2>/dev/null || exit 1
fi

case "${test_group}" in
  regression|slow)
    case "${compiler}" in
      gcc)
        case "${build}" in
          debug)
            export FFLAGS="-O0 -g -fcheck=all -fbacktrace -fbounds-check -ffpe-trap=invalid,zero,overflow -finit-real=nan -finit-integer=nan -Wall -Wunderflow -Wextra -Wconversion -Wuninitialized -Warray-bounds -Wshadow -Wno-unused-function -Wno-unused-variable -Wno-unused-parameter -Wno-unused-label -Wno-unused-but-set-variable"
            export OMP_NUM_THREADS=1
            ;;
          opt)
            export FFLAGS="-O1 -fopenmp -funroll-loops -finline-functions -ftree-vectorize -fstack-protector-strong"
            export OMP_NUM_THREADS=2
            ;;
          *)
            echo "Unknown build type: ${build}"
            return 1 2>/dev/null || exit 1
            ;;
        esac
        ;;
      intel|intel-classic)
        case "${build}" in
          debug)
            export FFLAGS="-O0 -debug all -check all -warn all,nodec,interfaces -gen_interfaces -traceback -fpe0 -ftrapuv -init=snan,arrays -check bounds"
            export OMP_NUM_THREADS=1
            ;;
          opt)
            export FFLAGS="-O -qopenmp -unroll -finline-functions -inline-forceinline -ipo -ip"
            export OMP_NUM_THREADS=2
            ;;
          *)
            echo "Unknown build type: ${build}"
            return 1 2>/dev/null || exit 1
            ;;
        esac
        ;;
      lfortran)
        case "${build}" in
          debug)
            export FFLAGS=""
            export OMP_NUM_THREADS=1
            ;;
          opt)
            export FFLAGS="--fast --openmp"
            export OMP_NUM_THREADS=2
            ;;
          *)
            echo "Unknown build type: ${build}"
            return 1 2>/dev/null || exit 1
            ;;
        esac
        ;;
      nvidia-hpc|flang)
        echo "${compiler} compiler not yet supported"
        return 1 2>/dev/null || exit 1
        ;;
      *)
        echo "Unknown compiler: ${compiler}"
        return 1 2>/dev/null || exit 1
        ;;
    esac
    ;;
  *)
    echo "Unknown test group: ${test_group}"
    return 1 2>/dev/null || exit 1
    ;;
esac
