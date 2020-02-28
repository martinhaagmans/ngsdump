import json
import glob
import pandas as pd

def json_to_excel(fn):
    with open(fn) as json_file:
        data = json.load(json_file)
        pd.DataFrame(data['molecularVariants']).to_excel(f'{fn}.xlsx')

for f in glob.glob('*.json'):
    json_to_excel(f)
