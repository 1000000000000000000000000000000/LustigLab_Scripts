import sys
import os
import requests

'''
	Downloads a single UniProt fasta file to a specified directory

	Example Command Line Use:

	python </path/to/download_single_uniprot_fasta_from_id.py> <UNIPROT ID> </path/to/desired/fasta/location>

'''

def download_single_uniprot_fasta_from_id(uniprot_id, out_path):
	url = f'https://www.uniprot.org/uniprot/{uniprot_id}.fasta'
	response = requests.get(url)
	if response.ok == True:
		print(f'Downloading {uniprot_id}.fasta.')
		fasta = response.text
		print(fasta)
		fasta_file = f'{uniprot_id}.fasta'
		fasta_file = os.path.join(out_path,fasta_file)
		with open(fasta_file,'w') as f:
    			f.write(response.text)
		return(response.text)
	else:
		print(f'Retrieving {uniprot_id} fasta failed!')
		return(None)

if __name__ == '__main__':

	uniprot_id = sys.argv[1]

	if len(sys.argv) == 3:
		out_path = sys.argv[2]
	else:
		out_path = os.getcwd()

	download_single_uniprot_fasta_from_id(uniprot_id, out_path)
