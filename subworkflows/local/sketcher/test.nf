#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { SKETCHER } from './main.nf' 

workflow test_sketcher {

    inputs = tuple( )

    SKETCHER ( inputs )
}
