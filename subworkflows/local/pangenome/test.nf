#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { PANGENOME } from './main.nf' 

workflow test_pangenome {

    inputs = tuple( )

    PANGENOME ( inputs )
}
