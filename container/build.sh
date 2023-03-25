#! /bin/sh
#
# build.sh
# Copyright (C) 2023 Tong LI <tongli.bioinfo@proton.me>
#
# Distributed under terms of the BSD-3 license.
#


docker build -t gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest .
singularity build /lustre/scratch126/cellgen/team283/imaging_sifs/gmm_decode.sif docker-daemon:gitlab-registry.internal.sanger.ac.uk/tl10/gmm-decoding:latest
