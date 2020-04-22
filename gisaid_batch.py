import os
from gisaid_record import GRec,GAlign,GVCFLine
from multiprocessing import Process,Manager
import numpy as np
from collections import defaultdict
import sys
import json
import re

class BatchProcess:
	""" methods to process batches of gisaid files producing more files of different types """

	def __init__(self,dir_path):
		self.dir_path = dir_path
		self.db_json_field_order = [
			'Accession ID', 'Virus name', 'Type', 'Passage details/history', 'Collection date', 'Comment', 'Location', 
			'Host', 'Additional location information', 'Gender', 'Patient age', 'Patient status', 'Specimen source', 
			'Additional host information', 'Outbreak', 'Last vaccinated', 'Treatment', 'Sequencing technology', 
			'Assembly method', 'Coverage', 'Originating lab', 'Address1', 'Sample ID given by the sample provider', 
			'Submitting lab', 'Address2', 'Sample ID given by the submitting laboratory', 'Authors', 'Submitter', 
			'Submission Date', 'Address3'
		]

	def _convert_html_to_json(self):
		""" convert all record html files to json files in the given directory """
		html_files = [name for name in os.listdir(self.dir_path) if name.endswith('.html')]
		for h in html_files:
			rec_id = h.rstrip('.html')
			html = ""
			with open(self.dir_path + '/' + h,'r') as f:
				html = f.read()
			gr = GRec(rec_id,html)
			gr.write_data_to_json_outfile()

	def _calculate_fasta_length(self,json_fasta):
		""" return the length of the fasta from a json record """
		return len(json_fasta.rstrip('\n').split('\n')[1])

	def _calculate_n_perc(self,json_fasta):
		""" return the percentage of Ns in the fasta from a json record """
		seq = json_fasta.rstrip('\n').split('\n')[1].upper()
		n_count = sum([1 for s in seq if s == 'N'])
		return round(100*n_count/len(seq),2)		

	def _clean_age_data_for_db(self,age):
		""" convert the age value into a float or null data type for db entry """
		age = age.replace(' ','')
		age = age.replace('\'','')
		if age.isnumeric():
			return age
		decimal_match = re.match(r'^\d+\.\d+$',age)
		if decimal_match is not None:
			return age
		decade_match1 = re.match(r'(\d)0s',age)
		decade_match2 = re.match(r'(\d)0-(\d)0',age)
		if (decade_match1 is not None) or (decade_match2 is not None):
			match = decade_match1
			if (decade_match2 is not None):
				match = decade_match2
			age = match.group(1) + '5'
			return age
		year_match = re.search(r'^(\d+)(\-)?((year)|\,).+',age) 
		if year_match is not None:
			age = year_match.group(1)
			return age
		month_match = re.match(r'^(\d+)months',age)
		if month_match is not None:
			age = round(int(month_match.group(1))/12,3)
			return age
		day_match =  re.match(r'^(\d+)days',age)			
		if day_match is not None:
			age = round(int(day_match.group(1))/365.25,3)
			return age
		return '\\N'

	def _clean_coverage_data_for_db(self,coverage):
		""" convert the coverage value into a float or null data type for db entry """
		coverage = coverage.replace(',','')
		coverage = coverage.replace('>','')
		coverage = coverage.replace('<','')
		x_match = re.match(r'([0-9\.]+)x',coverage)
		if x_match is not None:
			return x_match.group(1)
		x_match2 = re.search(r'x([0-9\.]+)',coverage)
		if x_match2 is not None:
			return x_match2.group(1)
		coverage = coverage.replace(' ','')
		coverage = coverage.rstrip('x')
		return '\\N'

	def _convert_json_file_to_db_line(self,json_filepath):
		""" given a json_filepath, make a tab-delimited line ready for entry into a DB """	
		with open(json_filepath) as f:
			data = json.load(f)
		seq_length = self._calculate_fasta_length(data['Fasta'])
		n_perc = self._calculate_n_perc(data['Fasta'])
		data['Patient age'] = self._clean_age_data_for_db(data['Patient age'])
		old = data['Coverage']
		data['Coverage'] = self._clean_coverage_data_for_db(data['Coverage'])
		line = "\t".join([str(data[key]).replace('\n',' ') for key in self.db_json_field_order]) + '\t' + str(seq_length) + '\t' + str(n_perc)
		return line

	def make_tab_delimited_db_file_from_json(self,outfile_path):
		""" from all the json files in a directory, make a single tab-delimited file ready for entry into a db table """
		json_files = [self.dir_path + '/' + name for name in os.listdir(self.dir_path) if name.endswith('.json')]
		with open(outfile_path,'w') as o:
			for jf in json_files:
				o.write(self._convert_json_file_to_db_line(jf) + '\n')

	def _make_vcf_files_from_html_file_list(self,html_file_list):
		""" helper method to make vcf files from a list of html file paths """
		for h in html_file_list:
			rec_id = h.rstrip('.html')
			html = ""
			with open(self.dir_path + '/' + h,'r') as f:
				html = f.read()
			gr = GRec(rec_id,html)
			ga = GAlign(gr.get_SeqRecord())
			ga.process_record()

	def threaded_make_vcf_files_from_html(self,num_threads,replace=False):
		""" produce a vcf file for each json record in the given directory, using multithreads """
		if not replace:
			vcf_files = set([name for name in os.listdir(self.dir_path) if name.endswith('.vcf')])
			html_files = [name for name in os.listdir(self.dir_path) if name.endswith('.html') and name not in vcf_files]
		else:
			html_files = [name for name in os.listdir(self.dir_path) if name.endswith('.html')]
		splits = np.array_split(html_files,num_threads)
		processes = []
		for i in range(len(splits)):
			p = Process( target=self._make_vcf_files_from_html_file_list, args=([splits[i]]) )
			processes.append(p)
			p.start()
	
	def get_vcf_line_record_dict_from_vcf_files(self,include_Ns=False):
		""" iterate through the vcf files in this directory and return a dict of vcf line : list of all records containing the line """
		vcf_lines = defaultdict(list)
		vcf_files = [name for name in os.listdir(self.dir_path) if name.endswith('.vcf')]
		for v in vcf_files:
			grec_id = v.rstrip('.vcf')
			with open(self.dir_path + '/' + v) as f:
				for line in f:
					if line.startswith('#') or line.startswith('CHROM'):
						continue
					if not include_Ns:
						alt = line.rstrip('\n').split('\t')[4]
						if 'N' in alt:
							continue
					vcf_lines[line.rstrip('\n')].append(grec_id)
		return vcf_lines
			
	def make_tab_delimited_db_files_from_vcf_lines(self,vcf_line_dict,var_outfile_path,record_outfile_path):
		""" given a list of variant lines with record info, make a tab-delimited file suitable for DB import """
		v = open(var_outfile_path,'w') 
		r = open(record_outfile_path,'w') 
		for var,records in vcf_line_dict.items():
			if ('\t0\t' in var):
				continue
			gv = GisVar(GVCFLine.from_line(var),records).get_var_type_child()
			v_line = "\t".join([gv.make_genomic_variant_name(),str(gv.genome_start),gv.get_type(),str(gv.is_synonymous(gv.make_protein_variant_name())),str(gv.get_protein_name()),str(gv.make_protein_variant_name())]) + '\n'
			v.write(v_line)
			for rec_id in records:
				r_line = rec_id + '\t' + gv.make_genomic_variant_name() + '\n'
				r.write(r_line)
		v.close()
		r.close()
				

if __name__ == '__main__':
	records_dir_path = sys.argv[1]
	bp = BatchProcess(records_dir_path)
	bp.make_tab_delimited_db_file_from_json('tmp')
	#bp.threaded_make_vcf_files_from_html(60,replace=False)
	from gisaid_variant import GisVar
	variants = bp.get_vcf_line_record_dict_from_vcf_files()
	bp.make_tab_delimited_db_files_from_vcf_lines(variants,'tmp2','tmp3')
	#for var,records in variants.items():
	#	if (len(records) < 2) or ('\t0\t' in var) or ('AAA\t' in var):
	#		continue
	#	gv = GisVar(GVCFLine.from_line(var),records).get_var_type_child()		
	#	print("\t".join([gv.get_type(),str(gv.is_synonymous(gv.make_protein_variant_name())),str(gv.get_protein_name()),str(gv.make_genomic_variant_name()),str(gv.make_protein_variant_name()),str(len(records)),records[0]]))
