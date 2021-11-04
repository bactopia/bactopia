#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { TBPROFILER } from './main.nf' 

workflow test_tbprofiler {

    inputs = tuple( )

    TBPROFILER ( inputs )
}
