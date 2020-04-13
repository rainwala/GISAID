import os
from gisaid_record import GRec,GAlign,GVCFLine
from multiprocessing import Process,Manager
import numpy as np
from collections import defaultdict
import sys

class BatchProcess:
	""" methods to process batches of gisaid files producing more files of different types """

	def __init__(self,dir_path):
		self.dir_path = dir_path
		pass

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

	def threaded_make_vcf_files_from_html(self,num_threads):
		""" produce a vcf file for each json record in the given directory, using multithreads """
		vcf_files = set([name for name in os.listdir(self.dir_path) if name.endswith('.vcf')])
		html_files = [name for name in os.listdir(self.dir_path) if name.endswith('.html') and name not in vcf_files]
		splits = np.array_split(html_files,num_threads)
		processes = []
		for i in range(len(splits)):
			p = Process( target=self._make_vcf_files_from_html_file_list, args=([splits[i]]) )
			processes.append(p)
			p.start()
	
	def get_variants_with_records_from_vcf_files(self,include_Ns=False):
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
			

if __name__ == '__main__':
	records_dir_path = sys.argv[1]
	bp = BatchProcess(records_dir_path)
	#bp.threaded_make_vcf_files_from_html(60)
	from gisaid_variant import GisVar
	variants = bp.get_variants_with_records_from_vcf_files()
	for var,records in variants.items():
		if (len(records) < 2) or ('\t0\t' in var) or ('AAA\t' in var):
			continue
		gv = GisVar(GVCFLine.from_line(var),records)		
		print(len(records))
		print(gv._get_genomic_variant_name())
