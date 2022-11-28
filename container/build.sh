#! /bin/sh
#
# build.sh
# Copyright (C) 2022 Tong LI <tongli.bioinfo@protonmail.com>
#
# Distributed under terms of the BSD-3 license.
#


docker build -t gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest .
singularity build /lustre/scratch117/cellgen/team283/imaging_sifs/gmm_decode.sif docker-daemon://gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest
