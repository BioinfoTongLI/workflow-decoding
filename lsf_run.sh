#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

MOUNT_POINT='/lustre/scratch126/cellgen/team283/NXF_WORK/'

TMP_NF_WORK="$MOUNT_POINT/${USER}_decoding_work"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK nextflow -trace nextflow.executor run BioinfoTongLI/workflow-decoding \
	-r develop \
	-params-file $1 \
	-with-report \
	-profile lsf \
	-resume
	#-entry $2 \
