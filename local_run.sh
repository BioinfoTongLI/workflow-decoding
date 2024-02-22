#! /bin/sh
#
# run.sh
# Copyright (C) 2021 Tong LI <tongli.bioinfo@gmail.com>
#
# Distributed under terms of the BSD-3 license.
#

MOUNT_POINT='/tmp/work/'

DATE_WITH_TIME=`date "+%Y%m%d%H%M"`
TMP_NF_WORK="$MOUNT_POINT/"

NXF_OPTS='-Dleveldb.mmap=false' NXF_VER=23.04.1 NXF_WORK=$TMP_NF_WORK nextflow -trace nextflow.executor run BioinfoTongLI/workflow-decoding \
	-r develop \
	-params-file $1 \
	-profile local \
	-with-report \
	-resume
