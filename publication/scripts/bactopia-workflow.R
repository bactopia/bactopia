library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

color_schema = list(
    default='#f7f7f7',
    optional='#fc8d59',
    user_dataset='#67a9cf',
    custom_dataset='#67a9cf'
)

workflow <- paste0("digraph workflow {
    # a 'graph' statement
    graph [
        overlap = false,
        layout = dot,
        nodesep = 2,
        ranksep = 5
    ]
    
    # several 'node' statements
    node [
        shape = box,
        fontname = Helvetica,
        fontsize = 144,
        height = 6,
        margin = 0.5,
        penwidth = 12,
        style = filled,
        fillcolor = '",color_schema$default,"'
    ]

    edge [arrowsize = 6, penwidth = 12]
    
    # Original FASTQ Input
    ofq1[label = 'Gather FASTQs']
    
    ofq2[label = 'Validate FASTQs']
    ofq4[label = 'Original Summary']
    ofq3[label = 'Genome Size']
    qc[label = 'Quality Control']
    ofq1->ofq2;
    ofq2->ofq3;
    ofq2->ofq4;
    ofq3->qc;
    
    # QC'd FASTQ Input
    a0[label = 'De novo Assembly']
    qc1[label = 'QC Summary']
    kmer[label = 'Count 31-mers']
    anno[label = 'Genome Annotation']
    maps[label = 'Reference Mapping', shape = 'oval', fillcolor='",color_schema$optional,"']
    snp[label = 'Call Variants', shape = 'oval', fillcolor='",color_schema$optional,"']
    amr[label = 'Antimicrobial Resistance']
    mlst[label = 'Sequence Type', shape = 'oval', fillcolor='",color_schema$optional,"']
    m0[label = 'Minmer Sketch']
    local[label = 'Local Assembly', shape = 'oval', fillcolor='",color_schema$optional,"']
    m1[label = 'Minmer Query', shape = 'oval', fillcolor='",color_schema$optional,"']
    blast[label = 'BLAST Alignment', shape = 'oval', fillcolor='",color_schema$optional,"']
    a0->mlst[arrowhead=dot];
    qc->mlst[arrowhead=dot];
    qc->{qc1,kmer,a0,m0}
    qc->{snp,maps,local}[arrowhead=dot];
    
    # Assembly Input
    a0->anno;
    anno->amr;
    empty1[style='invis'];
    empty2[style='invis'];
    #empty3[style='invis']
    amr->empty1[style='invis'];
    empty1->empty2[style='invis'];
    # empty2->empty3[style='invis'];
    
    aqc[label = 'Assembly Quality Control'];
    # Minmer Sketch Input
    
    
    # Variants Inputs
    node [shape = circle,
    fixedsize = true,
    width = 11]
    snp_auto[label = 'Nearest\nReference', fillcolor='",color_schema$custom_dataset,"']
    sketch[label = 'Public\nSketches', fillcolor='",color_schema$custom_dataset,"']
    m0->{m1,snp_auto}[arrowhead=dot];
    
    # Ariba Databases
    ariba[label = 'Ariba\nDatasets', fillcolor='",color_schema$custom_dataset,"']
    
    # Annotation Input
    protein_set[label = 'Custom\nProtein Set', fillcolor='",color_schema$custom_dataset,"']
    mlst_schema[label = 'MLST\nSchema', fillcolor='",color_schema$custom_dataset,"']
    # User Inputs
    node [height = 12,width=12, shape=diamond, style = filled, fillcolor = '",color_schema$user_dataset,"']
    snp_user[label = 'Reference']
    genes[label='Gene'];
    proteins[label='Protein'];
    primers[label='Primer']
    maps_user[label = 'Sequence']
    
    # Edges
    edge [arrowhead=box]
    ariba->local;
    snp_auto->snp;
    
    a0->blast;
    a0->aqc;
    snp_user->snp;
    primers->blast;
    proteins->blast;
    genes->blast;
    maps_user->maps;
    protein_set->anno;
    sketch->m1;
    mlst_schema->mlst;
    {rank=same;mlst}
    {rank=same;empty1;m1;snp}
    {rank=same;empty2;maps;blast;local}
    {rank=same;m0;a0;}
    {rank=sink;primers;proteins;genes;ariba;is_seqs; mlst_schema; maps_user; sketch;snp_user;}
}")

grViz(workflow)

grViz(workflow) %>%
  export_svg %>% charToRaw %>% rsvg_svg("bactopia-workflow-no-key.svg")
