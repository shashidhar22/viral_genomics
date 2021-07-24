#!/bin/bash

set -Eeuo pipefail

# Nextflow Version
NXF_VER=20.10.0

# Nextflow Configuration File
NXF_CONFIG=nextflow.config

# Workflow to Run (e.g. GitHub repository)
WORKFLOW_REPO=main.nf #shashidhar22/viral_genomics

# Queue to use for exeution
QUEUE=campus-new

# Load the Nextflow module (if running on rhino/gizmo)
ml Nextflow

# Load the Singularity module (if running on rhino/gizmo with Singularity)
ml Singularity
# Make sure that the singularity executables are in the PATH
export PATH=$SINGULARITYROOT/bin/:$PATH

# Run the workflow
NXF_VER=$NXF_VER \
nextflow \
    run \
    -c ${NXF_CONFIG} \
    ${WORKFLOW_REPO} \
    -entry removeHostAndAssemble \
    -w /fh/scratch/delete90/warren_h/ebv_enktl \
    -params-file main.yaml \
    -resume