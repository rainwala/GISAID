from gisaid_parser import GISAID
import sys
from datetime import date, timedelta
import time

from_date = date(2020,4,1)
to_date = date(2020,4,1)

g = GISAID(True)
print('object set up')
g.fill_date_form_field(g.form_field_id_dict['submission_date_from'],from_date)
time.sleep(2)
g.fill_date_form_field(g.form_field_id_dict['submission_date_to'],to_date)
g.process_records_for_all_pages()
