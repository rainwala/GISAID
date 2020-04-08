import os
from gisaid_record import GRec,GAlign

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

	def make_vcf_files_from_html(self):
		""" produce a vcf file for each json record in the given directory """
		html_files = [name for name in os.listdir(self.dir_path) if name.endswith('.html')]
		for h in html_files:
			rec_id = h.rstrip('.html')
			html = ""
			with open(h,'r') as f:
				html = f.read()
			gr = GRec(rec_id,html)
			ga = GAlign()
			ga.process_grec(gr)

if __name__ == '__main__':
	bp = BatchProcess('./')
	bp.make_vcf_files_from_html()
