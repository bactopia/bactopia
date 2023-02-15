#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { ASSEMBLER } from './main.nf' 

workflow test_assembler {

    inputs = tuple( )

    ASSEMBLER ( inputs )
}
