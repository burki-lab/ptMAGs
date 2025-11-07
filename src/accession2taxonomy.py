# a script that takes an input file with a list of 
# accession numbers and outputs a tsv with 
# column 1 = accession numbers
# column 2 = taxonomy string

from Bio import Entrez
import sys

def acc_to_taxo(input_file, output_file, database):
    '''
    Grab taxonomy from set of accession numbers
    '''
    # Provide your email address to NCBI
    Entrez.email = "mahwash.jamy@gmail.com"

    # Load your list of accession numbers from a file
    accession_numbers = []
    with open(input_file, "r") as file:
        accession_numbers = file.read().splitlines()

    # Initialize an empty list to store results
    results = []

    # Loop through the accession numbers and fetch the taxids
    for accession_number in accession_numbers:
        try:
            handle = Entrez.efetch(id=accession_number, db = database, retmode = "xml")
            record = Entrez.read(handle)

            # Get taxonomy string
            taxonomy = record[0].get("GBSeq_taxonomy")
            # Remove semicolons
            taxonomy = taxonomy.replace(';', '')

            # Taxonomy string doesnt include species name. So get organism too.
            organism = record[0].get("GBSeq_organism")

            # Put the two together (and remove the duplicate term)
            # Check if the last part of taxonomy matches the beginning of organism
            if taxonomy.endswith(organism.split(" ")[0]):
                # If there's a match, remove the overlapping part and join the strings
                taxonomy = " ".join([taxonomy, organism[len(organism.split(" ")[0]):]])
            else:
                # If there's no match, simply join the strings
                taxonomy = " ".join([taxonomy, organism])
            
            # Remove the double spaces
            taxonomy = " ".join(taxonomy.split())
            # Remove any periods
            taxonomy = taxonomy.replace('.', '')
            # Replace the spaces with underscores
            taxonomy = taxonomy.replace(' ', '_')
        except Exception as e:
            print(f"An error occurred for accession {accession_number}: {str(e)}")
            taxonomy = "N/A"
        results.append((accession_number, taxonomy))

    # Save the results to a tab-delimited file
    with open(output_file, "w") as outfile:
        for result in results:
            outfile.write("\t".join(result) + "\n")

    print("Done!")
