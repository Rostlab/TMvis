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

	publishDir "$output_path", mode:'copy', pattern: "*out.pdb", saveAs: {filename -> filename.replace("out.pdb","_ppm.PDB")}
	publishDir "$output_path", mode:'copy', pattern: "*datapar1"
	publishDir "$output_path", mode:'copy', pattern: "*datapar2"
	publishDir "$output_path", mode:'copy', pattern: "*datasub1", saveAs: {filename -> "${uniprot_id}.${filename}"}
   	output:
	/* save all PPM outputs
	output.pdb, datapar1, datapar2, datasub1
	*/
   		path "*"
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

        publishDir "$output_path", mode:'copy', pattern: "*out.pdb", saveAs: {filename -> filename.replace("out.pdb","_ppm.PDB")}
        publishDir "$output_path", mode:'copy', pattern: "*datapar1"
        publishDir "$output_path", mode:'copy', pattern: "*datapar2"
        publishDir "$output_path", mode:'copy', pattern: "*datasub1", saveAs: {filename -> "${uniprot_id}.${filename}"}

        output:
        /* save all PPM outputs
        output.pdb, datapar1, datapar2, datasub1
        */
                path "*"

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

fasta_datasubChannel_alpha = Channel.fromFilePairs("${input_dir_alpha}/**/*.{datasub1,fasta}", size: 2, checkIfExists: true, flat: true)
fasta_datasubChannel_beta = Channel.fromFilePairs("${input_dir_beta}/**/*.{datasub1,fasta}", size: 2, checkIfExists: true, flat: true)

// Annotate membrane location per residue in a FASTA format 
process membrane2fasta_alpha {
	errorStrategy 'ignore'
        input:
                tuple val(uniprot_id), path(datasub), path(fasta) from fasta_datasubChannel_alpha
        script:
                output_path = "${input_dir_alpha}/${uniprot_id}/${uniprot_id}_ppm3.fasta"
		"""
		#!python
		from Bio import SeqIO
		# need to save fasta files in format:
		# id, sequence, prediction per protein
		# file named with uniprot ID
		# read datasub1
		with open("${datasub}", 'r') as file:
			# e.g. 1(  22-  27), 2(  28-  33)
			lines_lst = file.readlines()
			datasub1 = lines_lst[0]

			# if datasub1 is empty -> protein is not a TM -> this returns an error
			# -> nextflow skips this particular protein and reports that the script was not completed
			# this way, we can identify by which proteins TM prediction has not worked 
			# natively in the nextflow report 

		ppm_membr_residues_str = datasub1.split(";")[3]
		ppm_membr_residues_str_lst = ppm_membr_residues_str.split(",")
		# clean list of numbers of residues which are inside of membrane
		# "m" - inside
		# "." - not transmembrane
		ppm_membr_residues_lst = []
		for membr_interval in ppm_membr_residues_str_lst:
			 # delete index and brackets
			interval_string = membr_interval.strip().strip(" )").split("(")[1]
			interval_lst = [int(i_str.strip("( )")) for i_str in interval_string.split("-")]
			# fortran uses 1 indexing, python 0
			membr_res = list(range(interval_lst[0]-1, interval_lst[1]-1))
			ppm_membr_residues_lst.extend(membr_res)
		print(ppm_membr_residues_lst)
		seq_obj = list(SeqIO.parse("${fasta}", "fasta"))[0]
		seq_len = len(seq_obj.seq)
		uniprotid = seq_obj.id.split('|')[1]

		middle = seq_len // 2
		prot_seq = seq_obj.seq[0:middle]
		membr_seq = seq_obj.seq[middle:seq_len]

		with open("${output_path}", "w") as fasta_file:
			writer = SeqIO.FastaIO.FastaWriter(fasta_file)
			# save protein entry as fasta
			ppm_membr_seq = ""
			for i, res in enumerate(membr_seq):
				if i in ppm_membr_residues_lst:
					ppm_membr_seq = ppm_membr_seq + "m"
				else:
					ppm_membr_seq = ppm_membr_seq + "."
			seq_obj.seq = prot_seq + ppm_membr_seq
			writer.write_record(seq_obj)
			print("File saved ", "${output_path}")
		"""

}

process membrane2fasta_beta {
	errorStrategy 'ignore'
	input:
		tuple val(uniprot_id), path(datasub), path(fasta) from fasta_datasubChannel_beta
	script:
		output_path = "${input_dir_beta}/${uniprot_id}/${uniprot_id}_ppm3.fasta"
		"""
		#!python
		from Bio import SeqIO
		# need to save fasta files in format:
		# id, sequence, prediction per protein
		# file named with uniprot ID
		# read datasub1
		with open("${datasub}", 'r') as file:
			# e.g. 1(  22-  27), 2(  28-  33)
			datasub1 = list(file.readlines())[0] # can through an error if not a TM, see description in membrane2fasta_alpha
		ppm_membr_residues_str = datasub1.split(";")[3]
		ppm_membr_residues_str_lst = ppm_membr_residues_str.split(",")
		# clean list of numbers of residues which are inside of membrane
		# "m" - inside
		# "." - not transmembrane
		ppm_membr_residues_lst = []
		for membr_interval in ppm_membr_residues_str_lst:
			 # delete index and brackets
			interval_string = membr_interval.strip().strip(" )").split("(")[1]
			interval_lst = [int(i_str.strip("( )")) for i_str in interval_string.split("-")]
			# fortran uses 1 indexing, python 0
			membr_res = list(range(interval_lst[0]-1, interval_lst[1]-1))
			ppm_membr_residues_lst.extend(membr_res)
		print(ppm_membr_residues_lst)
		seq_obj = list(SeqIO.parse("${fasta}", "fasta"))[0]
		seq_len = len(seq_obj.seq)
		uniprotid = seq_obj.id.split('|')[1]

		middle = seq_len // 2
		prot_seq = seq_obj.seq[0:middle]
		membr_seq = seq_obj.seq[middle:seq_len]
		
		with open("${output_path}", "w") as fasta_file:
			writer = SeqIO.FastaIO.FastaWriter(fasta_file)
			# save protein entry as fasta
			ppm_membr_seq = ""
			for i, res in enumerate(membr_seq):
				if i in ppm_membr_residues_lst:
					ppm_membr_seq = ppm_membr_seq + "m"
				else:
					ppm_membr_seq = ppm_membr_seq + "."
			seq_obj.seq = prot_seq + ppm_membr_seq
			writer.write_record(seq_obj)
			print("File saved ", "${output_path}")
		"""
}
