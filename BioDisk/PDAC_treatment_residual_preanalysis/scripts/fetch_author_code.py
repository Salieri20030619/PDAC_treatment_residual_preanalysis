from pathlib import Path
import requests,json,hashlib,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
PROXY={'http':'http://127.0.0.1:7897','https':'http://127.0.0.1:7897'}
records=[]
for repo in ['karthikj89/humanpdac','whwanglab/PDAC']:
    u='https://api.github.com/repos/'+repo
    r=requests.get(u,proxies=PROXY,timeout=(20,90));r.raise_for_status();branch=r.json()['default_branch']
    r=requests.get(u+'/git/trees/'+branch+'?recursive=1',proxies=PROXY,timeout=(20,90));r.raise_for_status();data=r.json()
    name=repo.replace('/','_');(ROOT/'00_docs'/f'{name}_tree.json').write_text(json.dumps(data,indent=2))
    print('REPO',repo,'COMMIT',data['sha'])
    candidates=[x for x in data['tree'] if x['type']=='blob' and x.get('size',0)<10_000_000 and any(t in x['path'].lower() for t in ['preprocess','metadata','readme','annotation','clinical','pdac_data','preprocessing'])]
    print('CANDIDATES',[(x['path'],x.get('size',0)) for x in candidates])
    for item in candidates[:15]:
        url=f'https://raw.githubusercontent.com/{repo}/{data["sha"]}/{item["path"]}'
        response=requests.get(url,proxies=PROXY,timeout=(20,90));response.raise_for_status()
        dest=ROOT/'00_docs'/('author_'+name+'_'+item['path'].replace('/','_'))
        dest.write_bytes(response.content);records.append({'url':url,'file':str(dest),'sha256':hashlib.sha256(response.content).hexdigest(),'bytes':len(response.content),'utc':datetime.datetime.now(datetime.timezone.utc).isoformat()})
(ROOT/'01_metadata/author_code_manifest.json').write_text(json.dumps(records,indent=2))
print('AUTHOR_CODE_FETCH_COMPLETE')
