#!/bin/bash

set -Eeuo pipefail

# Nextflow Version
NXF_VER=21.09.0-edge

# Nextflow Configuration File
NXF_CONFIG=nextflow.config

# Workflow to Run (e.g. GitHub repository)
WORKFLOW_REPO=main.nf #shashidhar22/viral_genomics
WORKFLOW_NAME=ebvAssembly
# Queue to use for exeution
QUEUE=campus-new

# Load the Nextflow module (if running on rhino/gizmo)
ml nextflow/21.09.0-edge

# Load the Singularity module (if running on rhino/gizmo with Singularity)
ml Singularity
# Make sure that the singularity executables are in the PATH
export PATH=$SINGULARITYROOT/bin/:$PATH

# Run the workflow
NXF_VER=$NXF_VER \
nextflow \
    run \
    ${WORKFLOW_REPO} \
    -c ${NXF_CONFIG} \
    -w /fh/scratch/delete90/warren_h/ebv_enktl \
    -entry ${WORKFLOW_NAME} \
    -with-tower \
    -params-file elshafa.yaml \
    -resume
    
