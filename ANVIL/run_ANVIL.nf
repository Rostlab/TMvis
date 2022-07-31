#!/usr/bin/env nextflow

work_dir = params.work_dir

input_dir = "${work_dir}/results/current/db/**/pLDDT90F1/**/*.pdb"
input_dir_alpha = "${work_dir}/results/current/db/alpha/pLDDT90F1"
input_dir_beta = "${work_dir}/results/current/db/beta/pLDDT90F1"

pdbChannel_alpha = Channel.fromPath("${input_dir_alpha}/**/*.pdb", checkIfExists: true)
pdbChannel_beta = Channel.fromPath("${input_dir_beta}/**/*.pdb", checkIfExists: true)
pdbChannel = Channel
		.fromPath("${input_dir}", checkIfExists: true)

process runANVIL_alpha {
	errorStrategy 'ignore'
	input:
		file(pdb) from pdbChannel_alpha

	publishDir "$output_path",  mode:'copy'
   	output:
   		file("${base}_ANVIL.pdb")
	script:
		base = pdb.baseName
                uniprot_id = base.split("-")[1]
                output_path = "${input_dir_alpha}/${uniprot_id}"
	
		"""
		python /TMvis/ANVIL/ANVIL.py -i $pdb -af -n /TMvis/ANVIL/naccess2.1.1/
		"""
}

process runANVIL_beta {
	errorStrategy 'ignore'
        input:
                file(pdb) from pdbChannel_beta

        publishDir "$output_path", mode:'copy'
        output:
                file("${base}_ANVIL.pdb")
        script:
                base = pdb.baseName
                uniprot_id = base.split("-")[1]
                output_path = "${input_dir_beta}/${uniprot_id}"

                """
                python /TMvis/ANVIL/ANVIL.py -i $pdb -af --beta -n /TMvis/ANVIL/naccess2.1.1/
                """
}
