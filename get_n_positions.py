""" gets the positions of all Ns in gisaid records, using the vcf and json files for each record, by looking for those vcf lines
that are exclusively Ns and represent a substitution of a reference sequence, and only in records above a minimum length cutoff """

import sys
import json
import os
import re
from collections import defaultdict

record_dir_path = sys.argv[1]
LEN_CUTOFF = 28000
REF_LEN = 29903
N_pattern = re.compile(r'^N+$')

N_location_counts = defaultdict(int)

json_files = [name for name in os.listdir(record_dir_path) if name.endswith('.json')]
vcf_files = [name for name in os.listdir(record_dir_path) if name.endswith('.vcf')]

def update_n_dict_from_vcf_file(vcf_filepath):
	with open(vcf_filepath) as v:
		for line in v:
			if line.startswith('#') or line.startswith('CHROM'):
				continue
			tabs = line.rstrip('\n').split('\t')
			pos = int(tabs[2])
			ref = tabs[3]
			alt = tabs[4]
			if len(alt) != len(ref):
				continue
			if N_pattern.match(alt) is None:
				continue
			for i in range(len(alt)):
				N_location_counts[pos] += 1
				pos += 1

for j in json_files:
	v = j.rstrip('.json') + '.vcf'
	if v not in vcf_files:
		continue
	j = record_dir_path + '/' + j
	v = record_dir_path + '/' + v
	## parse the json
	with open(j) as f:
		data = json.load(f)
	fasta_seq = data['Fasta'].split('\n')[1]
	if len(fasta_seq) < LEN_CUTOFF:
		continue
	## parse the vcf file for Ns	 
	update_n_dict_from_vcf_file(v)

for pos in range(1,REF_LEN+1):
	count = 0
	if pos in N_location_counts:
		count = N_location_counts[pos]
	print(pos,count)
