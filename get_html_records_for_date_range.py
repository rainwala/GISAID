from gisaid_parser import GISAID
import sys
from datetime import date, timedelta
import time

from_date = date(2020,4,12)
to_date = date(2020,4,12)

g = GISAID(headless=True,reverse_record_order=False)
form_field_values_dict = {
    'submission_date_from':from_date,
    'submission_date_to':to_date,  
  }
g.process_records_for_all_pages(form_field_values_dict)
