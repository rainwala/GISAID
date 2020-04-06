from gisaid_parser import GISAID
import sys
from datetime import date, timedelta
import time

from_date = date(2020,4,1)
to_date = date(2020,4,5)

g = GISAID(True)
print('object set up')
g.fill_date_form_field('submission_date_from',from_date)
time.sleep(2)
g.fill_date_form_field('submission_date_to',to_date)
time.sleep(2)
g.print_total_number_of_records_on_all_pages()
g.process_records_for_all_pages()
