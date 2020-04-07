from gisaid_record import GRec
import os

html_files = [name for name in os.listdir('./') if name.endswith('html')]

for h in html_files:
	if not h.endswith('.html'):
		continue
	rec_id = h.rstrip('.html')
	html = ""
	with open(h,'r') as f:
		html = f.read() 
	gr = GRec(rec_id,html)
	gr.write_data_to_json_outfile()
