# dssp_derived_vkabat.py
# Example Usage:
# python dssp_derived_vkabat.py --dir /path/to/dssp_csv_files --verbose

import pandas as pd
import numpy as np
import os
import argparse
import glob

def get_protein_list_from_csv_files(dir_path=os.getcwd()):
	# returns list of proteins from csv files
	proteins = []
	for file in os.listdir(dir_path):
		protein = file.split('_')
		if (file.split('.')[-1] == 'ipynb') or (file.split('.')[-1] == 'py'):
			pass
		else:
			if len(protein) == 1:
				pass
			elif len(protein) > 1:
				if protein[0][0] == '.':
					pass
				else:
					proteins.append(protein[0])
	print(list(np.unique(proteins)))
	return(list(np.unique(proteins)))

def make_dssp_vkabat_csv(protein_name, dir_path=os.getcwd()):
	# Input protein id
	# Output DF for the protein

	# compile list of the csv files needed to perform the calculations
	dssp_csv_files = glob.glob(protein_name + '*_dssp_acc.csv')

	base_df = pd.read_csv(dssp_csv_files[0]).drop(columns=['rawSS', 'Sec Struct', 'Sol Acc'])

	# Add each dssp column to right side of dataframe
	for i, file in enumerate(dssp_csv_files):
		# make series from each file containing the ss assignments
		# rename the column to the filename
		df = pd.read_csv(file)
		new_col_name = str(file.split('.csv')[0]) + '_SS'
		df.rename(columns={'Sec Struct':new_col_name}, inplace=True)
		#ss = df[str(file.split('.csv')[0])]
		ss = df[new_col_name]
		base_df = pd.concat([base_df, ss], axis="columns")

	# Calculate each term for DSSP-vkabat (ie. N, k, n1) and also calculate dssp-vkabat (k * N / n1)
	cols_to_count = base_df.columns[2:]

	# count the number of unique inputs and add to 'k' column
	base_df['dssp_k'] = base_df[cols_to_count].nunique(axis=1)

	# count the number of columns and add to 'N' column
	base_df['dssp_N'] = base_df[cols_to_count].count(axis=1)

	# count the most popular value and add to 'n1' column
	base_df['dssp_n1'] = base_df[cols_to_count].apply(lambda row: row.value_counts().iloc[0], axis=1)

	# calculate dssp_vkabat
	base_df['dssp_vkabat'] = base_df['dssp_k'] * base_df['dssp_N'] / base_df['dssp_n1']

	# make output csv
	out_csv_file_name = protein_name + '_dssp_vkabat.csv'
	base_df.to_csv(out_csv_file_name, index=False)

def make_acc_csv(protein_name, dir_path=os.getcwd()):
	# Input protein id
	# Output DF for the protein

	# compile list of the csv files needed to perform the calculations
	dssp_csv_files = glob.glob(protein_name + '*_dssp_acc.csv')

	base_df = pd.read_csv(dssp_csv_files[0]).drop(columns=['rawSS', 'Sec Struct', 'Sol Acc'])

	# Add each acc column to right side of dataframe
	for i, file in enumerate(dssp_csv_files):
		# make series from each file containing the ss assignments
		# rename the column to the filename
		df = pd.read_csv(file)
		new_col_name = str(file.split('.csv')[0])+'_ACC'
		df.rename(columns={'Sol Acc': new_col_name}, inplace=True)
		acc = df[new_col_name]
		base_df = pd.concat([base_df, acc], axis="columns")

	# Calculate each term for DSSP-vkabat (ie. N, k, n1) and also calculate dssp-vkabat (k * N / n1)
	cols_to_get_stats = base_df.columns[2:]

	# define a function to calculate mean, median, and standard deviation for each row
	def stats(row):
		mean = np.mean(row)
		median = np.median(row)
		std = np.std(row)
		return pd.Series([mean, median, std], index=['mean', 'median', 'std'])

	# apply the function to each row of the dataframe and add new columns for each statistic
	base_df[['mean', 'median', 'std']] = base_df[cols_to_get_stats].apply(stats, axis=1)
	base_df.rename(columns={'mean':'all_model_acc_mean', 'median':'all_model_acc_median', 'std':'all_model_acc_stdev'})

	# make output csv
	out_csv_file_name = protein_name + '_sol_acc.csv'
	base_df.to_csv(out_csv_file_name, index=False)

def main(dir_path=os.getcwd(), verbose=False):
	for protein in get_protein_list_from_csv_files(dir_path):
		if verbose == True:
			print(f'\nAbout to process files in: {dir_path}\n')
			print(f'Calculating dssp_vkabat for {protein}')
		make_dssp_vkabat_csv(protein, dir_path)
		if verbose == True:
			print(f'Calculating acc for {protein}')
		make_acc_csv(protein, dir_path)

if __name__ == '__main__':
	parser = argparse.ArgumentParser(
		description='This script is used to calculate DSSP-derived V_Kabat values and solvent accessibility summaries from *_dssp_acc.csv files.'
	)

	# add parser arguments
	parser.add_argument(
		'--dir',
		type=str,
		required=False,
		default=os.getcwd(),
		help='Path to the directory containing <protein name>*_dssp_acc.csv files. Defaults to the current working directory.'
	)

	parser.add_argument(
		'--verbose',
		action='store_true',
		help='Print additional progress messages.'
	)

	# parse the command line arguments
	args = parser.parse_args()

	# access the arguments using their names
	dir_path = os.path.abspath(args.dir)
	verbose = args.verbose

	# check that the input directory exists
	if not os.path.isdir(dir_path):
		print(f'Error: directory does not exist: {dir_path}')
		sys.exit(1)

	# change into the input directory
	# This is needed because make_dssp_vkabat_csv() and make_acc_csv()
	# use glob.glob() on the current working directory.
	os.chdir(dir_path)

	# run the script
	main(dir_path=os.getcwd(), verbose=verbose)
