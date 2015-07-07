#!/usr/bin/env bash

# Prepare the commands
PIPELINEDIR=/home/fuy2/data/piPipes
export BEDTOOLS=${PIPELINEDIR}/bin/bedtools_piPipes
export INSERT_BED_TO_BED2=${PIPELINEDIR}/bin/piPipes_insertBed_to_bed2
export PARAFLY=${PIPELINEDIR}/bin/ParaFly

export BED2LENDIS

# I added SolexaQA into repo/tools so this is no longer needed
# perl $PIPELINE_DIRECTORY/bin/SolexaQA_piPipes.pl
