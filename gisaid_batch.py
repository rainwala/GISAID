import os
from gisaid_record import GRec,GAlign
from multiprocessing import Process,Manager
import numpy as np

class BatchProcess:
	""" methods to process batches of gisaid files producing more files of different types """

	def __init__(self,dir_path):
		self.dir_path = dir_path
		pass

	def convert_html_to_json(self):
		""" convert all record html files to json files in the given directory """
		html_files = [name for name in os.listdir(self.dir_path) if name.endswith('.html')]
		for h in html_files:
			rec_id = h.rstrip('.html')
			html = ""
			with open(h,'r') as f:
				html = f.read()
			gr = GRec(rec_id,html)
			gr.write_data_to_json_outfile()

	def _make_vcf_files_from_html_file_list(self,html_file_list):
		""" helper method to make vcf files from a list of html file paths """
		for h in html_file_list:
			rec_id = h.rstrip('.html')
			html = ""
			with open(h,'r') as f:
				html = f.read()
			gr = GRec(rec_id,html)
			ga = GAlign()
			ga.process_grec(gr)

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
	

if __name__ == '__main__':
	bp = BatchProcess('./')
	bp.threaded_make_vcf_files_from_html(60)
