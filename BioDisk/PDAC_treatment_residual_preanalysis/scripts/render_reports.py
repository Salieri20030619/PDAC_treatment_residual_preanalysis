from pathlib import Path
import html,re,base64,json,hashlib,datetime
ROOT=Path('/mnt/BioDisk/wangmingyuan/PDAC_treatment_residual_preanalysis')
REPORT=ROOT/'08_reports'
def inline(s):
    s=html.escape(s)
    s=re.sub(r'`([^`]+)`',r'<code>\1</code>',s)
    s=re.sub(r'\*\*([^*]+)\*\*',r'<strong>\1</strong>',s)
    return re.sub(r'\[([^\]]+)\]\((https?://[^)]+)\)',r'<a href="\2">\1</a>',s)
def render(source):
    lines=source.splitlines();out=[];toc=[];i=0;heading_count=0
    while i<len(lines):
        line=lines[i].strip()
        if not line:i+=1;continue
        h=re.match(r'^(#{1,6}) (.*)',line)
        if h:
            level=len(h[1]);anchor=f'section-{heading_count}';heading_count+=1
            out.append(f'<h{level} id="{anchor}">{inline(h[2])}</h{level}>')
            if level==2:toc.append((anchor,h[2]))
            i+=1;continue
        im=re.fullmatch(r'!\[([^\]]*)\]\(([^)]+)\)',line)
        if im:
            p=(REPORT/im[2]).resolve();assert p.is_relative_to(ROOT)
            content=base64.b64encode(p.read_bytes()).decode()
            out.append(f'<figure><img alt="{html.escape(im[1])}" src="data:image/png;base64,{content}"></figure>');i+=1;continue
        if line.startswith('|'):
            rows=[]
            while i<len(lines) and lines[i].strip().startswith('|'):
                cells=[x.strip() for x in lines[i].strip().strip('|').split('|')]
                if not all(re.fullmatch('[-: ]+',x) for x in cells):rows.append(cells)
                i+=1
            assert all(len(row)==len(rows[0]) for row in rows),'Invalid Markdown table'
            out.append('<div class="tablewrap"><table><thead><tr>'+''.join('<th>'+inline(x)+'</th>' for x in rows[0])+'</tr></thead><tbody>')
            for row in rows[1:]:out.append('<tr>'+''.join('<td>'+inline(x)+'</td>' for x in row)+'</tr>')
            out.append('</tbody></table></div>');continue
        if line.startswith('- '):
            out.append('<ul>')
            while i<len(lines) and lines[i].strip().startswith('- '):out.append('<li>'+inline(lines[i].strip()[2:])+'</li>');i+=1
            out.append('</ul>');continue
        paragraph=[]
        while i<len(lines) and lines[i].strip() and not re.match(r'^(#|\||- |!\[)',lines[i].strip()):paragraph.append(lines[i].strip());i+=1
        if not paragraph:raise ValueError('Unsupported Markdown '+lines[i])
        out.append('<p>'+inline(' '.join(paragraph))+'</p>')
    nav='<nav aria-label="目录"><h2>目录</h2><ol>'+''.join(f'<li><a href="#{a}">{html.escape(t)}</a></li>' for a,t in toc)+'</ol></nav>'
    if len(toc)>3:out.insert(2,nav)
    return '\n'.join(out)
style='''body{max-width:1120px;margin:45px auto;padding:0 32px;color:#20252b;background:#fff;font:17px/1.85 "Microsoft YaHei","Noto Sans CJK SC",Arial,sans-serif}h1{font-size:32px;line-height:1.4;margin-bottom:14px}h2{font-size:24px;line-height:1.5;margin-top:45px;color:#15191d}h3{font-size:19px;line-height:1.5;margin-top:30px}p{margin:16px 0}nav{margin:32px 0;padding:4px 24px 18px;background:#f4f6f8}nav h2{margin-top:15px}nav ol{columns:2;column-gap:35px}nav li{break-inside:avoid;margin:6px 0}a{color:#245781;text-decoration:none}a:hover{text-decoration:underline}table{border-collapse:collapse;width:100%;font-size:14px;line-height:1.7;table-layout:auto}th,td{border:1px solid #cdd3d8;padding:12px 13px;text-align:left;vertical-align:top;overflow-wrap:anywhere}th{background:#f0f2f4;font-weight:600}thead{display:table-header-group}.tablewrap{overflow-x:auto;margin:24px 0}code{font-size:.87em;overflow-wrap:anywhere}figure{margin:30px 0}img{display:block;max-width:100%;height:auto}li{margin:9px 0}strong{color:#111}@media(max-width:700px){body{padding:0 18px;font-size:16px}nav ol{columns:1}h1{font-size:26px}table{font-size:12px}th,td{padding:8px}}@media print{@page{size:A4;margin:17mm}body{max-width:none;margin:0;padding:0;font-size:10pt;line-height:1.6}nav{display:none}h1{font-size:22pt}h2{font-size:15pt;break-after:avoid}h3{font-size:12pt;break-after:avoid}table{font-size:8pt}th,td{padding:5px}.tablewrap{overflow:visible}tr,figure{break-inside:avoid}p{orphans:3;widows:3}}'''
for stem in ['PDAC_client_research_proposal_zh','PDAC_preanalysis_zh','download_integrity_report']:
    md=(REPORT/(stem+'.md')).read_text(encoding='utf-8')
    assert '\ufffd' not in md
    body=render(md)
    title=md.splitlines()[0].lstrip('# ')
    output=f'<!doctype html><html lang="zh-CN"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>{html.escape(title)}</title><style>{style}</style></head><body>{body}</body></html>'
    (REPORT/(stem+'.html')).write_text(output,encoding='utf-8')
    print('RENDERED',stem,len(md),len(output))
print('CLIENT_REPORT_RENDER_COMPLETE')
