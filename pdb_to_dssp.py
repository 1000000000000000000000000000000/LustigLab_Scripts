#!/home/rpearson/.conda/envs/dssp/bin/python

# biopython is required
# For doc help visit:
#   https://www.biostars.org/p/427294/
#   https://biopython.org/docs/1.76/api/Bio.PDB.DSSP.html
#   https://biopython.org/wiki/The_Biopython_Structural_Bioinformatics_FAQ

import pandas as pd
from Bio.PDB import DSSP
from Bio.PDB import PDBList
from Bio.PDB.DSSP import dssp_dict_from_pdb_file
from Bio.PDB import PDBParser
from Bio.PDB.DSSP import make_dssp_dict
import os
import glob

def make_dssp_from_pdb_file(pdbFile):
    # Input pdbFile path to output a dssp dictionary
    dssp_tuple = dssp_dict_from_pdb_file(pdbFile)
    dssp_dict = dssp_tuple[0]
    return(dssp_dict)

def get_attributes_from_dssp_dict(input, mode='dict'):
    # Input variable should either be a dssp file or the output of the make_dssp_from_pdb_file with a pdb file input.
    if mode == 'file':
        # file should be some dssp file path, for example, "./5btrA.dssp"
        dssp = make_dssp_dict(file)
        parsed_dssp_return = list(dssp[0].items())
    elif mode == 'dict':
        dssp = input
        parsed_dssp_return = list(dssp.items())
    else:
        print("get get_attributes_from_dssp_dict function can accept mode='dict' or mode='file'.")

    res_pos_list = []
    res_list = []
    raw_ss_list = []
    sol_acc_list = []
    for return_tuple in parsed_dssp_return:
        raw_secondary_structure = return_tuple[1][1]
        res = return_tuple[1][0]
        res_position = return_tuple[0][1][1]
        acc = return_tuple[1][2]

        res_pos_list.append(res_position)
        res_list.append(res)
        raw_ss_list.append(raw_secondary_structure)
        sol_acc_list.append(acc)

    return(res_pos_list, res_list, raw_ss_list, sol_acc_list)

def simplify_raw_ss(raw_ss_list):
    # Input list of secondary structure assignments to output a
    # categorized or grouped secondary assignment list

    # DSSP ss conversion taken from table 6 in benjy's thesis where:
    # H, G, I --> Helix (H)
    # E and B --> Sheet (S)
    # T, S, '-' --> Other (O)
    
    simplified_ss_list = []
    for ss in raw_ss_list:
        if ss == 'H' or ss == 'G' or ss == 'I':
            simple_ss = 'H'
        if ss == 'E' or ss == 'B':
            simple_ss = 'S'
        if ss == 'T' or ss == 'S' or ss == '-':
            simple_ss = 'O'

        simplified_ss_list.append(simple_ss)
    return(simplified_ss_list)

def main(file):
    input = make_dssp_from_pdb_file(file)
    res_pos = get_attributes_from_dssp_dict(input, mode='dict')[0]
    residues = get_attributes_from_dssp_dict(input, mode='dict')[1]
    raw_ss = get_attributes_from_dssp_dict(input, mode='dict')[2]
    acc = get_attributes_from_dssp_dict(input, mode='dict')[3]
    simple_ss = simplify_raw_ss(get_attributes_from_dssp_dict(input, mode='dict')[2])

    data = {'Pos': res_pos, 'Res': residues, 'rawSS': raw_ss, 'Sec Struct': simple_ss, 'Sol Acc': acc}
    df = pd.DataFrame(data)

    if file.split('.')[-1] == 'dssp:':
        csv_file_name = file.split('.dssp')[0] + '_dssp_acc.csv'
    elif file.split('.')[-1] == 'pdb':
        csv_file_name = file.split('.pdb')[0] + '_dssp_acc.csv'
    else:
        print("Cound not make csv file. Make sure your input file is either pdb or dssp.")

    df.to_csv(csv_file_name, index=False)

if __name__ == '__main__':
    import sys

    if len(sys.argv) == 2:
        file = str(sys.argv[1])
        main(file)
    else:
        print("Please enter the dssp file path or the pdb file path.")
        print("Example:")
        print("python pdb_to_dssp_and_acc.py '/path/to/pdbfile.pdb'")
        print('\n---------------------------------------------------------------------------------------------\n')
        print(f'Proceeding with processing all pdb or dssp files in {os.getcwd()}')
        files = glob.glob('*.pdb')
        for file in files:
            print(f'Running {file}')
            main(file)
        print('Done.')
