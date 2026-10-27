from pathlib import Path
import subprocess,json,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
urls=[('ftp_https','https://ftp.ncbi.nlm.nih.gov/geo/series/GSE202nnn/GSE202051/suppl/GSE202051_totaldata-final-toshare.h5ad.gz'),('geo_download','https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE202051&format=file&file=GSE202051_totaldata-final-toshare.h5ad.gz')]
records=[]
for key,url in urls:
    cmd=['curl','--proxy','http://127.0.0.1:7897','-fL','--range','0-65535','--connect-timeout','15','--max-time','45','-D',str(ROOT/'logs'/f'{key}_probe.headers'),'-o',str(ROOT/'logs'/f'{key}_probe.bin'),url]
    p=subprocess.run(cmd,capture_output=True,text=True)
    rec={'url':url,'exit_code':p.returncode,'stderr':p.stderr[-1200:],'headers':(ROOT/'logs'/f'{key}_probe.headers').read_text(),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()};records.append(rec);print(json.dumps(rec),flush=True)
(ROOT/'logs/official_endpoint_probes.json').write_text(json.dumps(records,indent=2))
print('NETWORK_PROBE_COMPLETE')
