from tqdm import tqdm
from Bio import SeqIO
from Bio.PDB import *
import tarfile
import gzip
import glob
import os
import subprocess
import statistics


def get_secondary_structure(membr_seq):
    # count H, h, B, b, S and . and ignore the in/out topology (everything to lower case)
    membr_seq = membr_seq.lower()
    counts = {m: membr_seq.count(m) for m in set(membr_seq)}
    counts.pop('.')
    # determine type of protein: s, h or b or combination from a dictionary with counts for each per residue
    # assignment
    secondary_structure = ''.join(sorted(list(counts.keys())))
    return secondary_structure


def get_membrane_prots(proteome_path, convert2uniprot=False, prediction_line=True, check_alpha_beta_signal=False):
    records = list(SeqIO.parse(proteome_path, "fasta"))
    not_empty = 0
    alpha_beta_signal = {}
    membr_prot_IDs = []

    for seq_obj in records:
        # FASTA file contains ID, protein sequence and transmembrane residues prediction (3 lines)
        # len(protein sequence) = len(transmembrane prediction) = num(residues)
        # SeqIO parses everything from the second line after > as a sequence -> aa sequence and
        # transmembrane. prediction are captured together and needed to be split
        if prediction_line:
            seq_len = len(seq_obj.seq)
            middle = seq_len // 2
            prot_seq = seq_obj.seq[0:middle]
            membr_seq = seq_obj.seq[middle:seq_len]
            # effectively check if non-empty membrane prediction
            if not (membr_seq == len(membr_seq) * '.') and (len(prot_seq) == len(membr_seq)):
                not_empty += 1
                if check_alpha_beta_signal:
                    protein_type = get_secondary_structure(membr_seq)
                    if protein_type not in alpha_beta_signal.keys():
                        alpha_beta_signal[protein_type] = []
                        alpha_beta_signal[protein_type].append(seq_obj.id)
                    else:
                        alpha_beta_signal[protein_type].append(seq_obj.id)
                membr_prot_IDs.append(seq_obj.id)
        # normal FASTA file contains ID and protein sequence (2 lines)
        else:
            not_empty += 1
            membr_prot_IDs.append(seq_obj.id)
    print(not_empty, "/", len(records), "  predicted membrane proteins")

    if check_alpha_beta_signal:
        if convert2uniprot:
            alpha_beta_signal_uniprot = {}
            for key, value in alpha_beta_signal.items():
                # extract Uniprot ID: sp|Q5VST9|OBSCN_HUMAN -> Q5VST9
                uniprotIDs = [prot.split('|')[1] for prot in value]
                alpha_beta_signal_uniprot[key] = uniprotIDs
            return alpha_beta_signal_uniprot
        else:
            # return list of SeqIO IDs
            return alpha_beta_signal
    else:
        if convert2uniprot:
            # extract Uniprot ID: sp|Q5VST9|OBSCN_HUMAN -> Q5VST9
            uniprotIDs = [prot.split('|')[1] for prot in membr_prot_IDs]
            return uniprotIDs
        else:
            # return list of SeqIOs
            return membr_prot_IDs


def save_membrane_prots2fasta(proteome_path="./../data/current/human.txt", output_path="./../results/current/db/alpha/{}/",
                              exclude=None):
    # need to save fasta files in format:
    # id, sequence, prediction per protein
    # file named with uniprot ID
    records = list(SeqIO.parse(proteome_path, "fasta"))
    output_len = 0
    for seq_obj in records:
        uniprotid = seq_obj.id.split('|')[1]
        if (exclude is not None) and (uniprotid not in exclude):
            seq_len = len(seq_obj.seq)
            middle = seq_len // 2
            prot_seq = seq_obj.seq[0:middle]
            membr_seq = seq_obj.seq[middle:seq_len]
            if not (membr_seq == len(membr_seq) * '.') and (len(prot_seq) == len(membr_seq)):
                # create a directory for each protein
                os.makedirs(output_path.format(uniprotid), exist_ok=True)
                output_path_per_protein = output_path.format(uniprotid) + uniprotid + ".fasta"
                with open(output_path_per_protein, "w") as fasta_file:
                    # effectively check if non empty membrane prediction
                    writer = SeqIO.FastaIO.FastaWriter(fasta_file)
                    # save protein entry as fasta
                    seq_obj.seq = prot_seq + membr_seq
                    writer.write_record(seq_obj)
                    output_len += 1
                    #print("File saved ", output_path_per_protein)


#  Do not forget to Close the archive afterwards!
def read_AFdb(path):
    # pdb file names like 'AF-Q9NRQ5-F1-model_v2.pdb.gz'
    # cif files like 'AF-Q9NRR1-F1-model_v2.cif.gz'
    tar = tarfile.open(path, mode='r:')
    return tar


def getAFprots(path, convert2uniprot=False):
    tar = read_AFdb(path)
    afIDs = tar.getnames()
    print(len(afIDs), "  proteins in AFdb")
    tar.close()
    if convert2uniprot:
        uniprotIDs = [prot.split('-')[1] for prot in afIDs]
        return uniprotIDs
    else:
        return afIDs


def extractAFpdb(path, uniprot_ids, output_dir="./../results/current/db/alpha/{}"):
    tar = read_AFdb(path)
    af_file_names = tar.getnames()
    extract_af_file_names = []
    # get pathes to files which are predicted TMs
    for file_name in af_file_names:
        uniprot_id = file_name.split('-')[1]
        if uniprot_id in uniprot_ids:
            extract_af_file_names.append(file_name)
    # extract from archive and unzip
    for file_name in extract_af_file_names:
        file = tar.extractfile(file_name)
        uniprot_id = file_name.split('-')[1]
        with gzip.open(file, 'rt') as gzipfile:
            gunzipped = gzipfile.read()
            gunzipped_path = output_dir.format(uniprot_id) + '/' + file_name.strip('.gz')
            with open(gunzipped_path, 'wt') as file:
                file.write(gunzipped)
    print(len(extract_af_file_names), "AF files from ", path, "extracted to ", output_dir)
    tar.close()


# returns mean pLDDT of a protein structure
def getMeanpLDDT(pdb):
    parser = PDBParser()
    structure = parser.get_structure(pdb, pdb)
    pLDDTs = []
    assert len(structure) == 1, "More than 1 model in the pdb file"
    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    p = float(atom.get_bfactor())
                    pLDDTs.append(p)
                    break
    return statistics.mean(pLDDTs)


if __name__ == '__main__':
    root_dir = "./../"
    path_to_AFdb = os.path.join(root_dir, "data/current/UP000005640_9606_HUMAN_v2.tar")
    #path_to_AF_tm_subset = os.path.join(root_dir, "results/current/db")
    path_to_proteome = os.path.join(root_dir, "data/current/human.txt")
    # read uniprotIDs from AF archive
    af_proteins = getAFprots(path_to_AFdb, convert2uniprot=True)
    # read TMbed predictions: alpha, beta and signal
    alpha_beta_signal = get_membrane_prots(path_to_proteome,
                                           convert2uniprot=True,
                                           prediction_line=True,
                                           check_alpha_beta_signal=True)

    # alpha or helix ('h', 'hs')
    alpha_uniprotIDs = alpha_beta_signal['hs'] + alpha_beta_signal['h']
    beta_uniprotIDs = alpha_beta_signal['b']

    # check which TMbed proteins have AF structure
    notinAFdb = []
    for sec_str in ['hs', 'h', 'b']:
        predicted_membr_proteins = alpha_beta_signal[sec_str]
        if set(predicted_membr_proteins).issubset(set(af_proteins)):
            print("OK: there are AF 3D structures for all predicted {} membrane proteins".format(sec_str))
        else:
            # collect all proteins that are missing in AFdb no matter alpha, beta or mixed
            not_found_IDs = set(predicted_membr_proteins) - set(af_proteins)
            notinAFdb.extend(not_found_IDs)
            print(len(not_found_IDs), " predicted {} membrane proteins were not found in AFdb".format(sec_str))
            print(not_found_IDs)

    # separate alpha and beta to different directories
    db_dir = os.path.join(root_dir, "results/current/db/")
    db_alpha_dir = os.path.join(db_dir, "alpha/")
    db_beta_dir = os.path.join(db_dir, "beta/")
    os.makedirs(db_dir, exist_ok=True)
    os.makedirs(db_alpha_dir, exist_ok=True)
    os.makedirs(db_beta_dir, exist_ok=True)

    # save each protein in its own folder named by Uniprot ID
    db_alpha_ID_template = db_alpha_dir+"{}/"
    db_beta_ID_template = db_beta_dir+"{}/"

    # alpha helix ('h')
    save_membrane_prots2fasta(path_to_proteome, db_alpha_ID_template,
                              exclude=set(list(notinAFdb) +
                                          alpha_beta_signal['s'] +
                                          alpha_beta_signal['bhs'] +
                                          alpha_beta_signal['bh'] +
                                          alpha_beta_signal['b']))
    # beta barrel or strand ('b')
    save_membrane_prots2fasta(path_to_proteome, db_beta_ID_template,
                              exclude=set(list(notinAFdb) +
                                          alpha_beta_signal['s'] +
                                          alpha_beta_signal['bhs'] +
                                          alpha_beta_signal['bh'] +
                                          alpha_beta_signal['hs'] +
                                          alpha_beta_signal['h']))

    # extract all alpha and beta proteins from the AF archive
    extractAFpdb(path=path_to_AFdb, uniprot_ids=alpha_uniprotIDs,
                 output_dir=db_alpha_ID_template)
    extractAFpdb(path=path_to_AFdb, uniprot_ids=beta_uniprotIDs,
                 output_dir=db_beta_ID_template)

    # subset extraction
    # only extract those AF structures which
    # 1. have one file F1 (not longer than 1400 bp)
    # 2. have mean pLDDT > 90 (accurate structures)
    pLDDT_threshold = 90.0
    # counter
    pLDDT90F1 = []

    # create directories
    alpha_plddt_dir = os.path.join(db_alpha_dir, "pLDDT90F1")
    beta_plddt_dir = os.path.join(db_beta_dir, "pLDDT90F1")
    os.makedirs(alpha_plddt_dir, exist_ok=True)
    os.makedirs(beta_plddt_dir, exist_ok=True)
    # copy alpha and beta which pass F1 and pLDDT90 criteria to separate folders
    for sec_str in ['alpha', 'beta']:
        path_to_AF_tm_subset = db_dir + "{}/".format(sec_str)
        # extract IDs of one AF file proteins
        all_uniprot_ids_pdbs = [path.split('/')[-1].split('-')[1] for path in
                                glob.glob(os.path.join(path_to_AF_tm_subset, '**/*.pdb'), recursive=True)]
        counts_dict = {item: all_uniprot_ids_pdbs.count(item) for item in all_uniprot_ids_pdbs}
        one_file_af_membr_proteins = []
        for file_name in glob.glob(os.path.join(path_to_AF_tm_subset, '**/*.pdb'), recursive=True):
            uniprot_id = file_name.split('/')[-1].split('-')[1]
            if counts_dict[uniprot_id] == 1:
                one_file_af_membr_proteins.append(file_name)
        print("After F1 filtering # {} proteins: ".format(sec_str), len(one_file_af_membr_proteins))
        # copy relevant pdbs to a new folder
        for pdb in tqdm(one_file_af_membr_proteins):
            if getMeanpLDDT(pdb) > pLDDT_threshold:
                uniprot_id = pdb.split('/')[-1].split('-')[1]
                pLDDT90F1.append(uniprot_id)
                try:
                    subprocess.run(["mv", db_dir + '/' + sec_str + '/' + uniprot_id,
                                    db_dir + '/' + sec_str + '/' + 'pLDDT90F1' + "/" + uniprot_id])
                except:
                    subprocess.run(["mkdir", db_dir + '/' + sec_str + '/' + uniprot_id])
                    subprocess.run(["mv", db_dir + '/' + sec_str + '/' + uniprot_id,
                                    db_dir + '/' + sec_str + '/' + 'pLDDT90F1' + "/" + uniprot_id])
        print("mean pLDDT per {} protein > 90.0 found: ".format(sec_str), len(pLDDT90F1))
        pLDDT90F1 = []
