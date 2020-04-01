from splinter import Browser
import time
from datetime import date

class GISAID:
	""" methods to extract individual records from the GISAID website """
	def __init__(self,headless=True):
		self.url = 'https://www.epicov.org/epi3/frontend'
		self.username = 'rainwala'
		self.password = 'GuhRsX9o'
		self.form_field_id_dict = {
			'submission_date_from': 'ce_q83p7x_5e_input',
			'submission_date_to': 'ce_q83p7x_5f_input',
			'collection_date_from': 'ce_q83p7x_59_input',
			'collection_date_to': 'ce_q83p7x_5a_input',
		}
		self.browser = self._login(headless,1)
		self._get_to_covid_records_main_page(5)

	def _login(self,headless,sleep_time):
		""" login to the GISAID website """
		browser = Browser('firefox', headless=headless)
		browser.visit(self.url)
		time.sleep(sleep_time)
		browser.find_by_id('elogin').fill(self.username)
		browser.find_by_id('epassword').fill(self.password)
		browser.find_by_value('Login').first.click()
		return browser

	def _get_to_covid_records_main_page(self,wait_time):
		""" navigate to the covid-19 records page """
		self.browser.is_text_present('Browse',wait_time=wait_time)
		self.browser.find_by_value('Browse').click()

	def _get_record_div_list_from_current_page(self,wait_time):
		""" get the list of record divs from a given record page represented in the current browser state """
		self.browser.is_text_present('yui-dt-rec yui-dt-last yui-dt-even', wait_time=wait_time)
		divs = self.browser.find_by_css('td[headers="yui-dt0-th-e "] > div[class="yui-dt-liner"]')
		return divs

	def _write_record_html_from_div_list(self,div_list,sleep_time):
		""" for each record div in the div list, get the html, and return a dict of record ID : html """
		curr_div = 0
		record_html_dict = {}
		while curr_div < len(div_list):
			d = div_list[curr_div]
			try:
				d.click()
				time.sleep(sleep_time)
			## we need this exception here as it signifies that the frame html has been generated
			## (it is triggered when splinter tries to click on the next div behind the frame)
			except:
				time.sleep(sleep_time)
				## this next try catch block is for when the iframe hasn't loaded yet
				try:
					frame_elem = self.browser.find_by_css('iframe').first	
				except:
					continue
				## switch the context to the frame
				with self.browser.get_iframe(frame_elem) as frame:
					tds = frame.find_by_css('td')
					record_id = None
					for td in tds:
						if 'EPI_' in td.value:
							record_id = td.value
					if record_id is not None:
						print(record_id)
						html = frame.find_by_css('html').first.html
						## write the record via a GRec object
						gr = GRec(record_id,html)
						gr.write_html_to_outfile()
					curr_div += 1
					frame.find_by_css('img[src="/epi3/app_entities/entities/icons/24x24/navigate_left.png"]').first.click()

	def process_records_for_current_page(self):
		""" wrapper for the above methods to get parsed output file for each record from the current records page """
		div_list = self._get_record_div_list_from_current_page(5)
		self._write_record_html_from_div_list(div_list,1)

	def navigate_to_page(self,page_num):
		""" load the given page of records """
		try:
			page_elem = self.browser.find_by_css(f'a[page="{page_num}"]')
			if len(page_elem) == 0:
				return False
			page_elem.first.click()
			return True
		except:
			return False

	def fill_date_form_field(self,field_id,date):
		""" given a field_id and date object, fill out a date form field in the main records page 
		with a date in the format YYYY-MM-DD """
		form_elem = self.browser.find_by_css(f'input[id="{field_id}"]')
		if len(form_elem) == 0:
			return False
		datestring = date.strftime('%Y-%m-%d')
		form_elem.first.fill(datestring)
		return True
		
class GRec:
	""" data structures and methods to deal with individual records parsed from GISAID """
	def __init__(self,record_id,html):
		self.record_id = record_id
		self.html = html

	def write_html_to_outfile(self):
		""" write the record to the given outfile """
		outfile = self.record_id + ".html"
		with open(outfile,'w') as o:
			o.write(self.html)

if __name__ == "__main__":
	g = GISAID(False)
	#g.process_records_for_current_page()
	#more_pages = True
	#next_page = 2
	#while more_pages:	
	#	time.sleep(2)
	#	more_pages = g.navigate_to_page(next_page)
	#	next_page += 1
	date = date(2020,3,15)
	g.fill_date_form_field(g.form_field_id_dict['submission_date_from'],date)
	time.sleep(2)
	g.fill_date_form_field(g.form_field_id_dict['submission_date_to'],date)
