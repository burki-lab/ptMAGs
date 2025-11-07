# here we store functions to convert genbank files 
# to protein fasta files.
import os
from Bio import SeqIO

def refseq_gb_prot(input_dir, output_dir):
    '''
    Convert downloaded GenBank files to proteome fasta files.
    '''
    # Iterate over GenBank files in the input directory
    for genbank_file in os.listdir(input_dir):
        if genbank_file.endswith(".gb") or genbank_file.endswith(".gbk"):
            genbank_file_path = os.path.join(input_dir, genbank_file)
            
            # Open the GenBank file and iterate through the records
            for record in SeqIO.parse(genbank_file_path, "genbank"):
                accession = record.id
                # Remove the version number from the accession
                accession = accession.rsplit('.', 1)[0]
                proteins = []

                for feature in record.features:
                    # Check if the feature is a CDS (Coding Sequence) with a translation
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        protein_id = feature.qualifiers['protein_id'][0]
                        protein_seq = feature.qualifiers['translation'][0]
                        
                        # Include gene, accession, and protein_id in the header if "gene" is available
                        if 'gene' in feature.qualifiers:
                            gene = feature.qualifiers['gene'][0]
                            header = f">gene={gene};accession={accession};protein={protein_id}"
                            proteins.append(f"{header}\n{protein_seq}")

                # Determine the output file path in the output directory
                output_file_name = f"{accession}.fasta"
                output_fasta_file = os.path.join(output_dir, output_file_name)

                # Write the output FASTA file for the current genome
                if proteins:
                    with open(output_fasta_file, "w") as fasta_file:
                        fasta_file.write("\n".join(proteins))
                else:
                    # Create an empty FASTA file if no "gene" information is available
                    with open(output_fasta_file, "w") as fasta_file:
                        fasta_file.write('')


def mfannot_gb_prot(input_dir, output_dir, type):
    '''
    Convert MFannot GenBank files to proteome fasta files.
    '''
    # Iterate over GenBank files in the input directory
    for genbank_file in os.listdir(input_dir):
        if genbank_file.endswith(".gb") or genbank_file.endswith(".gbk"):
            genbank_file_path = os.path.join(input_dir, genbank_file)

            # Extract the accession from the file name
            accession = os.path.splitext(genbank_file)[0]
            
            # Initialise list of proteins
            proteins = []
            
            # Open the GenBank file and iterate through the records
            for record in SeqIO.parse(genbank_file_path, "genbank"):
                contig = record.id
 
                for feature in record.features:
                    # Check if the feature is a CDS (Coding Sequence) with a translation
                    if feature.type == "CDS" and "translation" in feature.qualifiers:
                        gene = feature.qualifiers['gene'][0]
                        protein_seq = feature.qualifiers['translation'][0]
                        if type == "ref":
                            header = f">gene={gene};accession={accession};contig={contig}"
                        elif type == "mag":
                            header = f">gene={gene};mag={accession};contig={contig}"
                        proteins.append(f"{header}\n{protein_seq}")

                # Determine the output file path in the output directory
                output_file_name = f"{accession}.fasta"
                output_fasta_file = os.path.join(output_dir, output_file_name)

                # Write the output FASTA file for the current genome
                with open(output_fasta_file, "w") as fasta_file:
                    fasta_file.write("\n".join(proteins))
