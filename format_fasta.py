#!/usr/bin/python

''' 
Description of script:
This script reformats the defline of FASTA files to allow
IsUnstruct to work properly.

'''

import re
from pathlib import Path

RE_ID4_LEADING = re.compile(r'^>([A-Za-z0-9]{4})_')          # >XXXX_
RE_UNIPROT = re.compile(r'^>\w+\|([A-Za-z0-9]{6,10})\|')     # >db|ACCESSION|...
RE_FIRST4 = re.compile(r'^>[^A-Za-z0-9]*([A-Za-z0-9]{4})')   # first 4 alnum after '>'

def reformat_header_line(header: str, chain: str = "A") -> tuple[str, str]:
	"""
	Format FASTA defline: >xxxx_CHAIN| description
	- Prefer 4-char ID from:
	  (1) existing >XXXX_, (2) UniProt ACCESSION, (3) first 4 alphanumerics after '>'
	- Lowercase the 4-char ID
	- Use provided 'chain' (default 'A')
	- Description = text after first space in original header (if any)
	"""
	orig = header.rstrip('\n')

	if not orig.startswith('>'):
		# No header at all; synthesize
		new_header = f">abcd_{chain}|\n{orig}"
		return new_header + '\n', "header missing; synthesized default"

	# Extract a 4-char ID
	m = RE_ID4_LEADING.match(orig)
	if m:
		id4 = m.group(1).lower()
	else:
		m = RE_UNIPROT.match(orig)
		if m:
			id4 = m.group(1)[:4].lower()
		else:
			m = RE_FIRST4.match(orig)
			id4 = (m.group(1).lower() if m else "abcd")

	# Build description (everything after the first space, if present)
	parts = orig.split('>', maxsplit=1)
	desc = parts[1] if len(parts) > 1 else ""

	new = f">{id4}_{chain}|"
	if desc:
		new += f" {desc}"

	return new + '\n', f"normalized to >{id4}_{chain}|"

def process_fasta(path: Path, inplace: bool, chain: str = "A") -> None:
	with path.open('r', encoding='utf-8', errors='replace') as f:
		lines = f.readlines()

	if not lines:
		print(f"WARNING: {path.name} empty; skipped.")
		return

	new_header, msg = reformat_header_line(lines[0], chain=chain)
	lines[0] = new_header

	out_path = path if inplace else path.with_name(path.stem + '_modified.fasta')
	with out_path.open('w', encoding='utf-8') as f:
		f.writelines(lines)

	print(f"{path.name}: {msg} -> wrote to {out_path.name}")

def main(file: Path | None = None, dir: Path = Path("."), chain: str = "A", inplace: bool = False) -> None:
	"""
	If 'file' is provided, process that file.
	Otherwise process all *.fasta in 'dir' (non-recursive).
	Default behavior processes all FASTA files in current working directory and writes *_modified.fasta copies.
	"""
	if file is not None:
		if not file.exists() or file.suffix.lower() != ".fasta":
			print(f"ERROR: {file} not found or not a .fasta")
			return
		process_fasta(file, inplace=inplace, chain=chain)
		print("format_fasta has finished running.")
		return

	if not dir.is_dir():
		print(f"ERROR: {dir} is not a directory")
		return

	files = sorted(p for p in dir.iterdir() if p.is_file() and p.suffix.lower() == ".fasta")
	if not files:
		print(f"\n\nERROR: There are no FASTA files in {dir}!\n\n")
		return

	for p in files:
		process_fasta(p, inplace=inplace, chain=chain)

	print("format_fasta.py has finished running.")

if __name__ == "__main__":
	import argparse

	parser = argparse.ArgumentParser(
		description="Reformats FASTA defline to look like the following: '>abcd_A'."
		            "Default (no args): Processes all .fasta files in the current working directory and writes a new '*_modified.fasta' file."
	)
	parser.add_argument("-f", "--file", type=Path, help="Single FASTA file to process")
	parser.add_argument("-d", "--dir", type=Path, default=Path("."), help="Directory with FASTA files (non-recursive)")
	parser.add_argument("-c", "--chain", default="A", help="Chain letter to enforce when matching '_.|' (default: A)")
	parser.add_argument("--inplace", action="store_true", help="Edit files in place which overwrites the originals")

	args = parser.parse_args()
	# Pass parsed values to main()
	main(file=args.file, dir=args.dir, chain=args.chain, inplace=args.inplace)
