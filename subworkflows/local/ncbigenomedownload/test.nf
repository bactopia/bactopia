#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

include { NCBIGENOMEDOWNLOAD } from './main.nf' 

workflow test_ncbigenomedownload {
    NCBIGENOMEDOWNLOAD()
}
