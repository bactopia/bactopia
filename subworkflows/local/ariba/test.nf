#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ARIBA } from './main.nf' 

workflow test_ariba {

    inputs = tuple( )

    ARIBA ( inputs )
}
