#! /bin/sh
#
# build.sh
# Copyright (C) 2023 Tong LI <tongli.bioinfo@proton.me>
#
# Distributed under terms of the BSD-3 license.
#


docker build -t bioinfotongli/decoding:latest .
docker push bioinfotongli/decoding:latest
