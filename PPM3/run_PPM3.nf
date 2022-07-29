#!/usr/bin/env nextflow

work_dir = "$PWD"

input_dir = "${work_dir}/results/current/db/**/pLDDT90F1/**/*.pdb"
input_dir_alpha = "${work_dir}/results/current/db/alpha/pLDDT90F1"
input_dir_beta = "${work_dir}/results/current/db/beta/pLDDT90F1"

pdbChannel_alpha = Channel.fromPath("${input_dir_alpha}/**/*.pdb", checkIfExists: true)
pdbChannel_beta = Channel.fromPath("${input_dir_beta}/**/*.pdb", checkIfExists: true)
pdbChannel = Channel
		.fromPath("${input_dir}", checkIfExists: true)

process runPPM3_alpha {
	input:
		file(pdb) from pdbChannel_alpha

	publishDir "$output_path", overwrite:true, mode:'copy'
   	output:
   		file("${base}out.pdb")
	script:
		base = pdb.baseName
                uniprot_id = base.split("-")[1]
                output_path = "${input_dir_alpha}/${uniprot_id}"
		//type of the membrane
		membrane="ERm"
		//use heteroatoms?
		hetat="0"
		//topology of N-terminus of first subunit included in the corresponding input PDB file
		topology=" in"
	
		"""
		cp /res.lib .
		cp /immers .
		
		# fit one membrane type per protein
		touch 1membrane.inp
		# number of membranes
		echo '1' >> 1membrane.inp
		# create a config file for PPM3
		echo " ${hetat} ${membrane} ${topology} ${pdb}" >> 1membrane.inp	
		./immers<1membrane.inp
		"""
}

process runPPM3_beta {
        input:
                file(pdb) from pdbChannel_beta

        publishDir "$output_path", overwrite:true, mode:'copy'
        output:
              file("${base}out.pdb")
        script:
                base = pdb.baseName
                uniprot_id = base.split("-")[1]
                output_path = "${input_dir_beta}/${uniprot_id}"
                //type of the membrane
                membrane="ERm"
                //use heteroatoms?
                hetat="0"
                //topology of N-terminus of first subunit included in the corresponding input PDB file
                topology=" in"

                """
                cp /res.lib .
                cp /immers .
                # fit one membrane type per protein
                touch 1membrane.inp
                # number of membranes
                echo '1' >> 1membrane.inp
                # create a config file for PPM3
                echo " ${hetat} ${membrane} ${topology} ${pdb}" >> 1membrane.inp
                ./immers<1membrane.inp
                """
}
