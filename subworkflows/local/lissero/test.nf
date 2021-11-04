#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { LISSERO } from './main.nf' 

workflow test_lissero {

    inputs = tuple( )

    LISSERO ( inputs )
}
