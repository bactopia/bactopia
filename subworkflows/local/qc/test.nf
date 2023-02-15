#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { QC } from './main.nf' 

workflow test_qc {

    inputs = tuple( )

    QC ( inputs )
}
