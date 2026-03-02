#Sample Database from Thera-SAbDab (Therapeutic Structural Antibody Database)
#By Oxford Protein Informatics Group

import pathlib
import requests

URL = "http://opig.stats.ox.ac.uk/webapps/newsabdab/static/downloads/TheraSAbDab_SeqStruc_OnlineDownload.csv"
OUT = pathlib.Path("data")
OUT.mkdir(exist_ok=True)
OUTFILE = OUT / "TheraSAbDab.csv"
 
r = requests.get(URL, timeout=60)
r.raise_for_status()
OUTFILE.write_bytes(r.content)

print("Saved:", OUTFILE.resolve())
