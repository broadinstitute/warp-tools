import pandas as pd
import requests
def get_new_gs_loc(nimp_aliquot_id):
	response = requests.get('https://assets.nemoarchive.org/api/alt-id/nhash/sample/' + nimp_aliquot_id)
	response.raise_for_status()  # Raises an error for 4xx/5xx responses
	data = response.json()
	gsLoc = []
	gsSize = []
	gsDate = []     # Parse the response as JSON
	for fi in data['files']:
		response2 = requests.get(fi)
		response2.raise_for_status()  # Raises an error for 4xx/5xx responses
		data2 = response2.json()
		gs_files = [f for f in data2['manifest_file_urls'] if f.get('url', '').startswith('gs://')]
		gsLoc.append(gs_files[0]['url'])
		gsSize.append(gs_files[0]['size'])
		gsDate.append(data2['last_modified'])
	gsTable=pd.DataFrame({'size':gsSize,'date':gsDate,'url':gsLoc})
	return gsTable