library(DiagrammeR)
library(DiagrammeRsvg)
library(rsvg)

color_schema = list(
    default='#f7f7f7',
    optional='#fc8d59',
    user_dataset='#67a9cf',
    custom_dataset='#67a9cf'
)

workflow <- paste0("digraph {
    graph [overlap = false,
           layout = dot, 
           nodesep = 0.25,
           ranksep = 0.1,
           labelloc = t,
           fontname = Helvetica]
    splines=false
    rankdir=LR
    size='8,1'
    ratio='compress'
    node [shape=plaintext, fontname = Helvetica, style = filled, fontsize=11]
    subgraph cluster1 {
        label = 'Key' ;
        color = black ;

        subgraph cluster_analysis {
            label = 'Analysis Condition'
            default_analysis [label='   Always Enabled   ',
                              shape = box,
                              fillcolor = '",color_schema$default,"'] ;
                          
            optional_analysis [label='   Dataset Enabled   ',
                               shape = 'oval',
                               fillcolor='",color_schema$optional,"'] ;
            default_analysis -> optional_analysis [style=invis] ;
        }
          
        subgraph cluster_types {
            label = 'Workflow Path'
            node[width=0.8]
            a [style=invis];
            b [style=invis] ;
            c [style=invis] ;
            d [style=invis] ;
            e [style=invis] ;
            f [style=invis] ;
            a -> b [label='Always Executed'] ;
            c -> d [label='Requires Additional Datasets', arrowhead=dot];
            e -> f [label='Dataset Input', arrowhead=box];
        }
        
        subgraph cluster_input {        
            label = 'Optional Dataset Inputs'
            custom_dataset [group=left, label='Public',
                            shape = circle, 
                            fillcolor='",color_schema$custom_dataset,"'] ;
                           
            user_dataset [group=left, label = 'User',
                          shape = diamond, 
                          height = 1,
                          width = 1,
                          fillcolor='",color_schema$user_dataset,"'];
    
            custom_dataset -> user_dataset [style=invis] ;
        }
        
    }
}")

grViz(workflow)

grViz(workflow) %>%
    export_svg %>% charToRaw %>% rsvg_svg("bactopia-key.svg")