#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

MOUNT_POINT='/tmp/work/'

DATE_WITH_TIME=`date "+%Y%m%d%H%M"`
TRACE_FILE="$MOUNT_POINT/registration_trace_${DATE_WITH_TIME}.tsv"
TMP_NF_WORK="$MOUNT_POINT/"

NXF_OPTS='-Dleveldb.mmap=false' NXF_WORK=$TMP_NF_WORK nextflow -trace nextflow.executor run /lustre/scratch117/cellgen/team283/tl10/gmm-decoding/main.nf \
	-with-trace $TRACE_FILE \
	-params-file $1 \
	-entry $2 \
	-profile local \
	-resume
