from splinter import Browser
from gisaid_record import GRec
import time
from datetime import date
from bs4 import BeautifulSoup
import json
import re
import os
import math

class GISAID:
	""" methods to extract individual records from the GISAID website """
	def __init__(self,headless=True,reverse_record_order=False,record_dir_path="/home/aawan/SHARED/COVID-19/GISAID/records/20200412"):
		self.url = 'https://www.epicov.org/epi3/frontend'
		self.username = 'rainwala'
		self.password = 'GuhRsX9o'
		self.RECORDS_PER_PAGE = 27 # the maximum number of records on a gisaid records page 
		self.form_field_position_dict = {
			'submission_date_from': 2,
			'submission_date_to': 3,
			'collection_date_from': 0,
			'collection_date_to': 1,
		}
		self.reverse_record_order = reverse_record_order
		self.record_dir_path = record_dir_path
		self.downloaded_record_set = self._get_downloaded_record_set()
		self.browser = self._login(headless,4)
		self._get_to_covid_records_main_page(4)

	def _get_downloaded_record_set(self):
		""" make a set of records that have already been downloaded in the record dir """
		return set([name.rstrip('.html') for name in os.listdir(self.record_dir_path) if (name.endswith('.html')) and (name.startswith('EPI'))])

	def _login(self,headless,wait_time):
		""" login to the GISAID website """
		browser = Browser('firefox', headless=headless)
		browser.visit(self.url)
		browser.is_text_present('elogin',wait_time=wait_time)
		browser.find_by_id('elogin').fill(self.username)
		browser.find_by_id('epassword').fill(self.password)
		browser.find_by_value('Login').first.click()
		print('logged in')
		return browser

	def _get_to_covid_records_main_page(self,wait_time):
		""" navigate to the covid-19 records page """
		self.browser.is_text_present('Browse',wait_time=wait_time)
		self.browser.find_by_value('Browse').click()
		print('on records page')

	def _get_page_count_from_record_count(self,record_count):
		""" get the total number of pages for the current record result set, based on the record count """
		return math.ceil(record_count/self.RECORDS_PER_PAGE)
	
	def _get_record_div_list_from_current_page(self,wait_time):
		""" get the list of record divs from a given record page represented in the current browser state """
		self.browser.is_text_present('yui-dt-rec yui-dt-last yui-dt-even', wait_time=wait_time)
		divs = self.browser.find_by_css('td[headers="yui-dt0-th-f "] > div[class="yui-dt-liner"]')
		new_divs = [d for d in divs if d.value not in self.downloaded_record_set]
		return new_divs

	def _check_for_captcha(self,frame,sleep_time):
		""" within a frame context, test to see if the frame is the capctcha pop-up, 
		and if so, prompt for user input to enter the captcha """
		input_elems = frame.find_by_css('input[maxlength="64"]')
		if len(input_elems) == 0:
			return False
		## output the captcha image url 
		image_url = self._get_captcha_image_url(frame)
		with open('image_urls','a') as i:
			i.write(image_url + "\n")
		print(f'The image url is: {image_url}')
		print("I await your command :)")
		captcha = input()
		print(f'You entered: {captcha}')
		input_elems.first.fill(captcha)
		time.sleep(sleep_time)
		ok_elem = frame.find_by_value('OK').first
		ok_elem.click()
		time.sleep(sleep_time)
		return True

	def _get_captcha_image_url(self,frame):
		""" inside a captcha frame, download the captcha image """
		image_elems = frame.find_by_css('div[class="sys-form-fi-captcha sys-event-hook"] > img')
		if len(image_elems) == 0:
			return None
		src = image_elems.first['src']
		return src

	def _write_record_from_div_list(self,div_list,sleep_time):
		""" for each record div in the div list, get the html, convert to json, and write the html and json """
		curr_div = 0
		record_html_dict = {}
		while curr_div < len(div_list):
			d = div_list[curr_div]
			while True:
				try:
					time.sleep(sleep_time)
					print(f'trying to click record {d.value}')
					d.click()
				## we need this exception here as it signifies that the frame html has been generated
				## (it is triggered when splinter tries to click on the next div behind the frame)
				except:
					time.sleep(sleep_time)
					## this next try catch block is for when the iframe hasn't loaded yet
					try:
						print('looking for record iframe')
						frame_elem = self.browser.find_by_css('iframe').first	
					except:
						continue
					## switch the context to the frame
					with self.browser.get_iframe(frame_elem) as frame:
						## check for the captcha first
						if self._check_for_captcha(frame,2):
							continue
						tds = frame.find_by_css('td')
						record_id = None
						for td in tds:
							if td.value.startswith('EPI_'):
								record_id = td.value
						if record_id is not None:
							print(record_id)
							html = frame.find_by_css('html').first.html
							## write the record via a GRec object
							gr = GRec(record_id,html)
							gr.write_html_to_outfile()
							gr.write_data_to_json_outfile()
						curr_div += 1
						frame.is_text_present('navigate_left',wait_time=10)
						frame.find_by_css('img[src="/epi3/app_entities/entities/icons/24x24/navigate_left.png"]').first.click()
						break
	
	def _reverse_records(self,wait_time=3):
		""" reverse the order of records in the main records page """
		self.browser.is_text_present('Accession ID',wait_time=wait_time)
		accession_header = self.browser.find_by_value('Accession ID').first
		accession_header.click()

	def _get_selected_page(self,wait_time=3):
		""" return the current, selected, page number """
		self.browser.is_text_present('yui-pg-current-page yui-pg-page',wait_time=wait_time)
		return int(self.browser.find_by_css('span[class="yui-pg-current-page yui-pg-page"]').first.value)
		
	def _process_records_for_current_page(self):
		""" wrapper for the above methods to get parsed output file for each record from the current records page """
		div_list = self._get_record_div_list_from_current_page(30)
		print(f'{len(div_list)} new records in this page')
		self._write_record_from_div_list(div_list,3)

	def _navigate_to_page(self,page_num,wait_time=10):
		""" load the given page of records """
		try:
			print(f'trying to go to page {page_num}')
			self.browser.is_text_present(f'page="{page_num}"',wait_time=wait_time)
			page_elem = self.browser.find_by_css(f'a[page="{page_num}"]')
			if len(page_elem) == 0:
				print(f'failed to go to page {page_num} as it is not clickable')
				return False
			page_elem.first.click()
			## make sure that we have navigated to the correct page
			navigated_page = self._get_selected_page()
			while navigated_page != page_num:
				navigated_page = self._get_selected_page()
			print(f'successfully gone to page {page_num}')
			return True
		except:
			return False

	def _fill_date_form_field(self,field_name,date,wait_time=15):
		""" given a field_id and date object, fill out a date form field in the main records page 
		with a date in the format YYYY-MM-DD """
		if field_name not in self.form_field_position_dict:
			return False
		field_position = self.form_field_position_dict[field_name]
		self.browser.is_text_present('sys-form-fi-date',wait_time=wait_time)
		form_elem = self.browser.find_by_css('div[class="sys-form-fi-date"] > input')[field_position]
		datestring = date.strftime('%Y-%m-%d')
		print('getting record count')
		record_count = self._get_total_number_of_records_on_all_pages()
		print('filling out date form')
		form_elem.fill(datestring)
		print('filled!')
		## click somewhere else to return focus to the records page
		print('clicking elsewhere')
		self.browser.is_text_present('Search',wait_time=wait_time)
		self.browser.find_by_value('Search').first.click()
		#form_elem.click()
		print('clicked!')
		## wait until the record count has changed before proceeding
		print('checking for updated record count')
		new_record_count = self._get_total_number_of_records_on_all_pages()
		while new_record_count == record_count:
			new_record_count = self._get_total_number_of_records_on_all_pages()
		return True

	def _filter_records(self,form_fields_values_dict,sleep_time=3):
		""" fill out the form to filter records based on the provided values """
		for field,value in form_fields_values_dict.items():
			print(f'filling out {field} with {value}')
			self._fill_date_form_field(field,value)
		if self.reverse_record_order: # sometimes we want to process records in reverse id order
			time.sleep(sleep_time)
			print('reversing records')
			self._reverse_records()
			time.sleep(sleep_time)

	def _get_total_number_of_records_on_all_pages(self,wait_time=15):
		""" get the total number of records found on all pages of the main records page """
		self.browser.is_text_present('white-space:nowrap',wait_time=wait_time)
		record_count = self.browser.find_by_css('span[style="white-space:nowrap;"]').last.value
		record_count = re.sub(r'^\s*Total:\s+(\S+)\s+.+',r'\1',record_count)
		record_count = record_count.replace(',','')
		return int(record_count)

	def process_records_for_page_range(self,form_fields_values_dict,start_page=1,end_page=False,sleep_time=3):
		""" get the html files for records for all pages for the given form field values, assuming a start at page 1 """
		## get the total number of records and the total number of pages
		self._filter_records(form_fields_values_dict,sleep_time=sleep_time)
		record_count = self._get_total_number_of_records_on_all_pages()
		page_count = self._get_page_count_from_record_count(record_count)
		if not end_page:
			end_page = page_count
		print(f'There are {record_count} records in {page_count} pages, and we are processing pages {start_page} to {end_page}')
		## process the first page
		print(f'page {start_page}')
		if start_page > 1:
			pages_per_display = 4
			current_last_page = 10
			while(current_last_page < start_page):
				self._navigate_to_page(current_last_page)
				current_last_page += pages_per_display
			self._navigate_to_page(start_page)
		self._process_records_for_current_page()
		## process the remaining pages
		for i in range(start_page + 1,end_page + 1):
			self._navigate_to_page(i)
			time.sleep(sleep_time)
			self._process_records_for_current_page()
		print('shutting down...')

		
if __name__ == "__main__":
	g = GISAID(True)
	date = date(2020,2,1)
	form_field_values_dict = {
		'submission_date_from':date,
		'submission_date_to':date,
	}
	g.process_records_for_all_pages(form_field_values_dict)
