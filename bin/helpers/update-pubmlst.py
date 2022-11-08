#!/usr/bin/env python
PROGRAM = "update-pubmlst.py"
DESCRIPTION = "Download latest XML from PubMLST"
VERSION = "2.1.1"

def get_pubmlst_xml():
    import requests
    r = requests.get("https://pubmlst.org/static/data/dbases.xml")
    return r.text


def parse_pubmlst_xml(text):
    """
    <species>Yersinia ruckeri
        <mlst>
            <database>
                <url>https://pubmlst.org/yruckeri/</url>
                <profiles>
                    <url>
                        https://rest.pubmlst.org/db/pubmlst_yruckeri_seqdef/schemes/1/profiles_csv
                    </url>
                </profiles>
                <loci>
                    <locus>glnA
                    <url>https://rest.pubmlst.org/db/pubmlst_yruckeri_seqdef/loci/glnA/alleles_fasta</url>
                    </locus>
                    <locus>gyrB
                    <url>https://rest.pubmlst.org/db/pubmlst_yruckeri_seqdef/loci/gyrB/alleles_fasta</url>
                    </locus>
                    <locus>dnaJ
                    <url>https://rest.pubmlst.org/db/pubmlst_yruckeri_seqdef/loci/dnaJ/alleles_fasta</url>
                    </locus>
                    <locus>thrA
                    <url>https://rest.pubmlst.org/db/pubmlst_yruckeri_seqdef/loci/thrA/alleles_fasta</url>
                    </locus>
                    <locus>hsp60
                    <url>https://rest.pubmlst.org/db/pubmlst_yruckeri_seqdef/loci/hsp60/alleles_fasta</url>
                    </locus>
                    <locus>recA
                    <url>https://rest.pubmlst.org/db/pubmlst_yruckeri_seqdef/loci/recA/alleles_fasta</url>
                    </locus>
                </loci>
            </database>
        </mlst>
    </species>    
    """
    import xml.etree.ElementTree as ET
    pubmlst = {}
    stage = "data"
    species = ""
    locus = ""
    root = ET.fromstring(text)
    for descendant in root.iter():
        if descendant.tag == "data":
            continue
        elif descendant.tag == "species":
            species = descendant.text.rstrip()
            locus = ""
            pubmlst[species] = {
                "db": "",
                "profiles": "",
                "loci": {}
            }
        elif descendant.tag == "mlst":
            continue
        elif descendant.tag == "database":
            stage = "db"
        elif descendant.tag == "profiles":
            stage = "profiles"
        elif descendant.tag == "loci":
            stage = "loci"
        elif descendant.tag == "locus":
            locus = descendant.text.rstrip()
        elif descendant.tag == "url":
            if stage == "db":
                pubmlst[species]["db"] = descendant.text.rstrip()
            elif stage == "profiles":
                pubmlst[species]["profiles"] = descendant.text.rstrip()
            elif stage == "loci":
                pubmlst[species]["loci"][locus] = descendant.text.rstrip()
    
    return pubmlst


if __name__ == '__main__':
    import argparse as ap
    import datetime
    import json
    import sys
    import textwrap

    parser = ap.ArgumentParser(
        prog=PROGRAM,
        conflict_handler='resolve',
        description=(
            f'{PROGRAM} (v{VERSION}) - {DESCRIPTION}'
        ),
        formatter_class=ap.RawDescriptionHelpFormatter
    )
    parser.add_argument('bactopia', metavar="STR", type=str,
                        help='Directory where Bactopia repository is stored.')
    parser.add_argument('--version', action='version',
                        version=f'{PROGRAM} {VERSION}')

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(0)

    args = parser.parse_args()
    pubmlst_xml = get_pubmlst_xml()
    pubmlst_dict = parse_pubmlst_xml(pubmlst_xml)

    with open(f'{args.bactopia}/data/pubmlst.json', 'wt') as fh_out:
        json.dump(pubmlst_dict, fh_out, sort_keys=True, indent=4)
    with open(f'{args.bactopia}/data/pubmlst-updated.txt', 'wt') as fh_out:
        fh_out.write(datetime.datetime.now().isoformat())
        fh_out.write("\n")
