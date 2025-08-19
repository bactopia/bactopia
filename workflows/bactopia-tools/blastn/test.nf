include { BLASTN } from './main.nf'
include { BLASTN as BLASTN_GZ } from './main.nf'
include { BLASTN as BLASTN_USE_GENES } from './main.nf'

workflow test_blastn {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTN(
        inputs,
        file(params.test_data['species']['portiera']['genome']['single_gene'], checkIfExists: true)
    )
}

workflow test_blastn_gz {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTN_GZ(
        inputs,
        file(params.test_data['species']['portiera']['genome']['single_gene_gz'], checkIfExists: true)
    )
}

workflow test_blastn_use_genes {

    inputs = tuple(
        [ id:"SRX1390622" ],
        file(params.test_data['species']['portiera']['genome']['blastdb'], checkIfExists: true)
    )

    BLASTN_USE_GENES(
        inputs,
        file(params.test_data['species']['portiera']['genome']['single_gene'], checkIfExists: true)
    )
}
