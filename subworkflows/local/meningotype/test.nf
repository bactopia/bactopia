#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { MENINGOTYPE } from './main.nf' 

workflow test_meningotype {

    inputs = tuple( )

    MENINGOTYPE ( inputs )
}
