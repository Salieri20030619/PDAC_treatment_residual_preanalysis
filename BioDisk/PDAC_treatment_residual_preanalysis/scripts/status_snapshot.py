from pathlib import Path
import json,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
out={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'downloads':{},'completed_samples':[],'refined_samples':[]}
for ds in ['GSE205013','GSE202051']:
    folder=ROOT/'02_raw'/ds/'segments';plan=json.loads((folder/'plan.json').read_text());parts=list(folder.glob('*.bin'))
    have=sum(p.stat().st_size for p in parts);out['downloads'][ds]={'bytes_downloaded':have,'bytes_total':plan['total_bytes'],'fraction':have/plan['total_bytes']}
out['completed_samples']=[p.name.split('_')[1].split('.')[0] for p in (ROOT/'03_processed').glob('GSE205013_P*.complete.json')]
out['refined_samples']=[p.name.split('_')[1] for p in (ROOT/'03_processed').glob('GSE205013_P*_cnv_refined.json')]
out['stages']={p.name:json.loads(p.read_text())['exit_code'] for p in (ROOT/'logs').glob('*.status.json')}
print(json.dumps(out,indent=2));(ROOT/'logs/latest_status_snapshot.json').write_text(json.dumps(out,indent=2))
print('STATUS_SNAPSHOT_COMPLETE')
