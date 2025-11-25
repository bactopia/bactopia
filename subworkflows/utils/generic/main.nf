//
// Generic functions to help save code duplication across multiple workflows
//
nextflow.preview.types = true

//
// Gather results from a channel
//
// Returns: Channel<Tuple<Map, Set<Path>>>
//
def gather(Channel<Tuple<Map, Path>> ch_results, String toolName) {
    return ch_results
        .collect{_meta, output -> output}
        .map{ output -> 
            tuple([id:toolName], output.toSet()) as Tuple<Map, Set<Path>>
        }
}


//
// Mix multiple channels and flatten file sets
// Tuple<Map, Set<Path>>  -->  Tuple<Map, Path>
//
// Returns: Channel<Tuple<Map, Path>>
//
def flattenPaths(List<Channel<Tuple<Map, Set<Path>>>> channels) {
    // Handle single channel case
    if (channels.size() == 1) {
        return channels[0].flatMap { meta, files -> 
            files.collect { file -> tuple(meta, file) }
        }
    }
    
    // Mix multiple channels and flatten
    def mixed = channels[0]
    channels[1..-1].each { ch -> mixed = mixed.mix(ch) }
    return mixed.flatMap { meta, files -> 
        files.collect { file -> tuple(meta, file) }
    }
}

//
// Adapt 4-element sample tuples to appropriate size based on data availability
// 
// Returns one of:
// - Channel<Tuple<Map, Set<Path>>>
// - Channel<Tuple<Map, Set<Path>, Set<Path>>>
// - Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>> 
//
def formatSamples(Channel<Tuple<Map, Set<Path>, Set<Path>, Set<Path>>> samples, Value<Integer> dataTypes) {
    if (dataTypes.value == 1) {
        return samples.map { meta, inputs, _extra, _extra2 -> 
            tuple(meta, inputs) as Tuple<Map, Set<Path>>
        }
    } else if (dataTypes.value == 2) {
        return samples.map { meta, inputs, extra, _extra2 -> 
            tuple(meta, inputs, extra) as Tuple<Map, Set<Path>, Set<Path>>
        }
    } else if (dataTypes.value == 3) {
        return samples
    } else {
        error("Invalid total data types: ${dataTypes}. Must be 1, 2, or 3")
    }
}
