"""Offline HybAs report from recorded outputs. No analysis tools are executed."""
import csv
import hashlib
import html
import io
import json
import math
from datetime import datetime, timezone
from pathlib import Path
from hybas_report_assets import CSS

SCHEMA = 'hybas-report-1'
RENDERER_VERSION = 'HybAs-9.0.1-AWARE'


def e(value):
    return html.escape(str(value), quote=True)


def number(value):
    try:
        value = float(str(value).replace(',', ''))
        return value if math.isfinite(value) else None
    except (TypeError, ValueError):
        return None


def show(value, unit=''):
    if value is None or value in ('', 'NA', 'N/A', 'unknown'):
        return 'Not recorded'
    if isinstance(value, bool):
        return 'Yes' if value else 'No'
    n = number(value)
    return (f'{n:,.0f}' if n is not None and n == int(n) else f'{n:,.2f}' if n is not None else e(value)) + unit


def table(rows, columns):
    return '<div class="table-wrap"><table><thead><tr>'+''.join('<th>'+e(c.replace('_',' '))+'</th>' for c in columns)+'</tr></thead><tbody>'+''.join('<tr>'+''.join('<td>'+show(row.get(c))+'</td>' for c in columns)+'</tr>' for row in rows)+'</tbody></table></div>'


def collect(root):
    """Capture bounded metadata, source checksums and annotation counts for one isolate."""
    root = Path(root)
    sources, warnings = {}, []
    def read(rel, required=False):
        path = root / rel
        if not path.is_file():
            if required:
                raise ValueError('Required saved output missing: '+rel)
            return None
        if path.stat().st_size > 20*1024*1024:
            warnings.append('Metadata too large to embed: '+rel)
            return None
        raw = path.read_bytes()
        sources[rel] = hashlib.sha256(raw).hexdigest()
        return raw.decode('utf-8-sig')
    summary = json.loads(read('metadata/run_summary.json', True))
    if not isinstance(summary, dict) or not summary.get('isolate') or not summary.get('workflow_version'):
        raise ValueError('Invalid HybAs run summary')
    tables = {}
    for name in ['kraken_decision','qc_assessment','skip_reasons','tool_versions','workflow_provenance','system_info','checksums']:
        text = read('metadata/'+name+'.tsv')
        tables[name] = list(csv.DictReader(io.StringIO(text), delimiter='\t')) if text else []
    benchmarks = []
    for path in sorted((root/'benchmarks').glob('*.tsv')):
        text = read(str(path.relative_to(root)))
        if text:
            benchmarks.extend(dict(rule=path.stem, **row) for row in csv.DictReader(io.StringIO(text), delimiter='\t'))
    tables['benchmarks'] = benchmarks
    annotation = {'status':'Not recorded','counts':{}}
    if summary.get('run_prokka') is False:
        annotation['status'] = 'Disabled in recorded configuration'
    else:
        path = root/'annotation'/(str(summary['isolate'])+'.gff')
        if path.is_file():
            counts, digest = {}, hashlib.sha256()
            with path.open('rb') as handle:
                in_fasta = False
                for raw in handle:
                    digest.update(raw)
                    if raw.startswith(b'##FASTA'):
                        in_fasta = True
                    if in_fasta or raw.startswith(b'#'):
                        continue
                    parts = raw.decode('utf-8', errors='replace').rstrip().split('\t')
                    if len(parts) == 9:
                        counts[parts[2]] = counts.get(parts[2],0)+1
            sources[str(path.relative_to(root))] = digest.hexdigest()
            annotation = {'status':'Recorded GFF feature counts','counts':counts}
    artifacts = []
    for label, rel in [('Final assembly','work/final/assembly.final.fasta'),('Prokka GFF','annotation/'+str(summary['isolate'])+'.gff'),('Protein FASTA','annotation/'+str(summary['isolate'])+'.faa'),('MultiQC','reports/multiqc/multiqc_report.html'),('QUAST table','reports/quast/report.tsv'),('fastp','reports/fastp/fastp.html')]:
        artifacts.append({'name':label,'path':rel,'present_at_capture':(root/rel).is_file() and (root/rel).stat().st_size>0})
    return dict(schema_version=SCHEMA, captured_at=datetime.now(timezone.utc).isoformat(), summary=summary,
                tables=tables, annotation=annotation, artifacts=artifacts, sources_sha256=sources, warnings=warnings)


def normalize(data):
    if data.get('schema_version') == SCHEMA and isinstance(data.get('summary'),dict):
        return data
    if data.get('isolate') and data.get('workflow_version'):
        return dict(schema_version=SCHEMA, summary=data, tables={}, annotation={}, artifacts=[], sources_sha256={},
                    warnings=['This edition uses a legacy run summary. Detailed platform cleaning actions, annotation counts and file inventories were not captured in that summary.'])
    raise ValueError('Unsupported HybAs report evidence')


def keyvalues(rows):
    result = {}
    for row in rows:
        values = list(row.values())
        if len(values)>=2:
            result[values[0]] = values[1]
    return result


def methods_text(bundle):
    s=bundle['summary']
    stages=s.get('executed_stages')
    stages=', '.join(stages) if isinstance(stages,list) else 'not recorded'
    return (f"Sample {s.get('isolate','not recorded')} was processed using {s.get('workflow_version','not recorded')} in {s.get('run_mode','not recorded')} mode. "
            f"The recorded final assembly source was {s.get('final_source_stage','not recorded')}. Executed polishing stages: {stages or 'none'}. "
            f"The requested polishing mode was {s.get('requested_polish_mode',s.get('polish_mode_requested','not recorded'))}. "
            f"Recorded AWARE decision mode: {s.get('kraken_decision_mode','not recorded')}; selected target: {s.get('kraken_selected_target_name','not recorded')} "
            f"(taxid {s.get('kraken_selected_target_taxid','not recorded')}). "
            f"BUSCO enabled: {s.get('run_busco','not recorded')}; lineage: {s.get('busco_lineage','not recorded')}. "
            f"Prokka enabled: {s.get('run_prokka','not recorded')}. "
            "This draft describes the saved run; review it alongside the recorded configuration and tool versions before publication.")


def render_report(data):
    b=normalize(data);s=b['summary'];ts=b.get('tables',{});decisions=keyvalues(ts.get('kraken_decision',[]))
    methods=methods_text(b)
    def stat(label,value,note=''):
        return '<div class="stat"><div class="stat-label">'+e(label)+'</div><div class="stat-value">'+show(value)+'</div><div class="stat-note">'+e(note)+'</div></div>'
    def section(ident,title,content):
        return '<section class="section" id="'+ident+'"><div class="section-head"><h2>'+title+'</h2></div>'+content+'</section>'
    def panel(title,content):
        return '<div class="panel"><h3>'+title+'</h3>'+content+'</div>'
    def detail(title,content):
        return '<details class="panel"><summary>'+title+'</summary>'+content+'</details>'
    def pairs(values):
        return '<dl class="kv">'+''.join('<div><dt>'+e(k)+'</dt><dd>'+show(v)+'</dd></div>' for k,v in values)+'</dl>'
    def bars(rows,unit):
        vals=[number(v) for _,v in rows];maximum=100 if unit=='%' else max([v for v in vals if v is not None]+[1])
        return '<div class="metric-bars">'+''.join('<div class="metric-row"><span>'+e(label)+'</span><div class="metric-track">'+('<span style="width:'+str(max(0,100*n/maximum))+'%"></span>' if n is not None else '')+'</div><strong>'+show(value,unit)+'</strong></div>' for (label,value),n in zip(rows,vals))+'</div>'
    notices=''.join('<p class="callout warn">'+e(w)+'</p>' for w in b.get('warnings',[]))
    overview=notices+'<div class="stats">'+stat('Assembly length',s.get('assembly_bp'),'bp · saved final assembly')+stat('Contigs',s.get('assembly_contig_count_fasta'))+stat('N50',s.get('quast_n50'),'bp · QUAST')+stat('BUSCO complete',s.get('busco_complete_percent') if s.get('run_busco') is not False else None,'Percent · disabled' if s.get('run_busco') is False else 'Percent · recorded lineage')+'</div>'
    overview+=panel('Run at a glance',pairs([('Input mode',s.get('run_mode')),('Final source stage',s.get('final_source_stage')),('GC (%)',s.get('quast_gc_percent')),('Analysis version',s.get('workflow_version')),('Report renderer',RENDERER_VERSION),('Circularity','Not recorded in this evidence bundle')]))
    long_prefix = "pacbio" if s.get("long_read_platform") == "pacbio" else "ont"
    long_title = "PacBio" if long_prefix == "pacbio" else "Oxford Nanopore"
    def read_metric(prefix, key, metric):
        if prefix == 'pacbio':
            return (s.get('pacbio') or {}).get(key, {}).get('bases' if metric == 'bp' else metric)
        return s.get(prefix+'_'+key+'_'+metric)
    read_panels=[]
    for title,prefix,stages in [('Illumina R1','illumina',[('Raw','raw_r1'),('After cleaning','clean_r1'),('After trimming','trim_r1')]),('Illumina R2','illumina',[('Raw','raw_r2'),('After cleaning','clean_r2'),('After trimming','trim_r2')]),(long_title,long_prefix,[('Raw','raw'),('After cleaning','cleaned')] if long_prefix=='pacbio' else [('Raw','raw'),('After cleaning','clean'),('After Filtlong','filt')])]:
        inactive=(prefix=='illumina' and s.get('run_mode') in {'ont_only','pacbio_only'}) or (prefix in {'ont','pacbio'} and s.get('run_mode')=='illumina_only')
        body='<p class="muted">Not used in this run.</p>' if inactive else bars([(label,read_metric(prefix,key,'reads')) for label,key in stages],' reads')+detail('Recorded base counts',table([dict(stage=label,bases=read_metric(prefix,key,'bp')) for label,key in stages],['stage','bases']))
        read_panels.append(panel(title,body))
    reads='<p class="muted">Each platform and mate is shown separately. Bars compare recorded read counts within that panel.</p><div class="grid-2">'+''.join(read_panels)+'</div>'
    aware=panel('Selected target',pairs([('Name',s.get('kraken_selected_target_name')),('Taxid',s.get('kraken_selected_target_taxid')),('Rank',s.get('kraken_selected_target_rank')),('Selection source',s.get('kraken_selected_target_source')),('Selection reason',s.get('kraken_target_selection_reason'))]))
    actions=[]
    for name,prefix in [('Illumina','illumina'),(long_title,long_prefix)]:
        inactive=(prefix=='illumina' and s.get('run_mode') in {'ont_only','pacbio_only'}) or (prefix in {'ont','pacbio'} and s.get('run_mode')=='illumina_only')
        actions.append(dict(platform=name,target_support_percent='Not used' if inactive else s.get('pacbio_selected_target_percent') if prefix=='pacbio' else s.get('kraken_'+prefix+'_selected_target_percent'),cleaning_action='Not used' if inactive else decisions.get(prefix+'_kraken_decision','Not recorded'),reason=decisions.get(prefix+'_kraken_decision_reason','Not recorded')))
    aware+=panel('Platform cleaning decisions',table(actions,['platform','target_support_percent','cleaning_action','reason']))
    aware+=detail('Confidence and cross-platform agreement',pairs([('Decision confidence',s.get('kraken_decision_confidence')),('AWARE state',s.get('kraken_awareness_state')),('Agreement score',s.get('kraken_platform_agreement_score')),('Agreement reason',s.get('kraken_platform_agreement_reason')),('Runner-up',s.get('kraken_runner_up_name'))])+'<p>Target selection, platform-specific cleaning actions and combined confidence are separate recorded decisions.</p>')
    executed=s.get('executed_stages',[]);executed=executed if isinstance(executed,list) else []
    skipped=s.get('skipped_stages',{});skipped=skipped if isinstance(skipped,dict) else {}
    stages=[dict(stage=s.get('assembler') or 'unicycler',status='Recorded assembly source',reason='')]+[dict(stage=x,status='Executed',reason='') for x in executed]+[dict(stage=x,status='Skipped',reason=reason) for x,reason in skipped.items()]
    assembly=panel('Recorded assembly route','<div class="stage-route">'+''.join('<span class="badge '+('ok' if row['status']=='Executed' else 'neutral')+'">'+e(row['stage'])+' · '+e(row['status'])+'</span>' for row in stages)+'</div>'+pairs([('Requested polishing',s.get('requested_polish_mode',s.get('polish_mode_requested'))),('Final source',s.get('final_source_stage')),('Final source SHA-256',s.get('final_source_sha256'))])+detail('Stage details',table(stages,['stage','status','reason'])))
    assembly+=panel('Final assembly metrics',table([{'metric':label,'value':s.get(key)} for label,key in [('Length (bp)','quast_total_length'),('Largest contig (bp)','quast_largest_contig'),('N50 (bp)','quast_n50'),('L50','quast_l50'),('GC (%)','quast_gc_percent')]],['metric','value'])+'<p class="small muted">Stage-by-stage quality measurements were not captured in this summary; no improvement is inferred from a stage having run.</p>')
    quality=panel('BUSCO',('<p>Disabled in the recorded configuration.</p>' if s.get('run_busco') is False else '<p>Lineage: '+show(s.get('busco_lineage'))+'</p>'+bars([(label,s.get('busco_'+key+'_percent')) for label,key in [('Single copy','single_copy'),('Duplicated','duplicated'),('Fragmented','fragmented'),('Missing','missing')]],'%')))
    quality+=panel('Coverage',('<p>PacBio mapped coverage is not recorded in this evidence bundle.</p>' if long_prefix=='pacbio' else pairs([('Estimated ONT coverage (×)',s.get('estimated_ont_coverage_x')),('Recorded ONT bases',s.get('coverage_ont_bp')),('Assembly bases used',s.get('coverage_assembly_bp'))]))+'<p class="small muted">This is the saved base-count estimate, not a per-position mapped-depth distribution.</p>')
    if ts.get('qc_assessment'):
        quality+=detail('Workflow QC checks',table(ts['qc_assessment'],list(ts['qc_assessment'][0])))
    annotation=b.get('annotation',{})
    annotation_html=panel('Recorded annotation', '<p>'+e(annotation.get('status','Not recorded'))+'</p>'+table([dict(feature=k,count=v) for k,v in annotation.get('counts',{}).items()],['feature','count'])+'<p>Use the final assembly, matching GFF and protein FASTA for downstream BaSIS and GeneScope runs.</p>')
    artifacts=b.get('artifacts',[])
    outputs=panel('Output inventory',table([dict(name=a['name'],path=a['path'],availability='Present at capture' if a['present_at_capture'] else 'Not available at capture') for a in artifacts],['name','path','availability'])+'<p class="small muted">Paths are relative to the isolate output folder. Open these files through Bactologix Results & downloads or the output folder. MultiQC remains a detailed tool report.</p>')
    versions=ts.get('tool_versions') or [dict(tool=k.removeprefix('tool_').removesuffix('_version'),version=v) for k,v in s.items() if k.startswith('tool_') and k.endswith('_version')]
    methodology=panel('Methods draft','<p>'+e(methods)+'</p><div class="download-list"><button class="button" id="download-methods">Download methods</button><button class="button" id="download-json">Download report evidence</button></div>')
    methodology+=detail('Tool versions',table(versions,list(versions[0]) if versions else ['tool','version']))
    for key,title in [('workflow_provenance','Recorded configuration'),('system_info','System information'),('checksums','Workflow output checksums'),('benchmarks','Recorded per-rule resource use')]:
        rows=ts.get(key,[])
        if rows:methodology+=detail(title,table(rows,list(rows[0])))
    methodology+=detail('Report evidence provenance',pairs([('Schema',b['schema_version']),('Captured',b.get('captured_at')),('Analysis timestamp',s.get('timestamp'))])+table([dict(path=k,sha256=v) for k,v in b.get('sources_sha256',{}).items()],['path','sha256']))
    blocks=[('overview','Assembly overview',overview),('reads','Reads & cleaning',reads),('aware','AWARE decisions',aware),('assembly','Assembly & polishing',assembly),('quality','Genome quality',quality),('annotation','Annotation & handoff',annotation_html),('outputs','Detailed outputs',outputs),('methods','Methods & provenance',methodology)]
    encoded=json.dumps(b,ensure_ascii=False).replace('<',r'\u003c').replace('>',r'\u003e').replace('&',r'\u0026')
    return '<!doctype html><html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1"><title>HybAs AWARE · '+e(s.get('isolate','Sample'))+'</title><style>'+CSS+'''\n.metric-bars{display:grid;gap:18px;margin:24px 0}.metric-row{display:grid;grid-template-columns:110px 1fr 120px;gap:12px;align-items:center;font-size:12px}.metric-row strong{text-align:right}.metric-track{height:12px;background:#e8edf1;border-radius:6px;overflow:hidden}.metric-track span{display:block;height:100%;background:#087c78;border-radius:6px}.stage-route{display:flex;gap:10px;flex-wrap:wrap;margin-bottom:22px}.kv dd{overflow-wrap:anywhere}.main{padding-top:25px}.hero{display:block}.section>p{margin-bottom:18px}@media(max-width:740px){.metric-row{grid-template-columns:90px 1fr;gap:8px}.metric-row strong{grid-column:2;text-align:left}.main{padding-top:17px}}\n</style></head><body><a class="skip-link" href="#overview">Skip to report</a><div class="app"><aside class="sidebar"><a class="brand" href="#overview"><span>HybAs<small>AWARE · BACTOMICS</small></span></a><p class="nav-label">Assembly report</p><nav aria-label="Report sections">'''+''.join('<a class="nav-link" href="#'+ident+'">'+e(title)+'</a>' for ident,title,_ in blocks)+'''</nav><div class="sidebar-foot">Offline report · saved run evidence</div></aside><main class="main"><header class="hero"><p class="eyebrow">BACTOMICS / HYBAS AWARE</p><h1>'''+e(s.get('isolate','Sample'))+'</h1><p class="hero-description">Assembly, cleaning decisions and genome quality</p><div class="hero-meta"><span>Mode: '+show(s.get('run_mode'))+'</span><span>Analysis: '+show(s.get('workflow_version'))+'</span><span>Renderer: '+RENDERER_VERSION+'</span></div></header>'+''.join(section(*block) for block in blocks)+'''<footer>Recorded measurements and decisions · unavailable evidence is not a passing result.</footer><script type="application/json" id="report-data">'''+encoded+'''</script><script>
const data=JSON.parse(document.getElementById('report-data').textContent);
function download(name,text,type){const url=URL.createObjectURL(new Blob([text],{type}));const a=document.createElement('a');a.href=url;a.download=name;a.click();setTimeout(()=>URL.revokeObjectURL(url),1000)}
document.getElementById('download-json').onclick=()=>download('hybas-report-evidence.json',JSON.stringify(data,null,2),'application/json');
document.getElementById('download-methods').onclick=()=>download('hybas-methods.txt',document.querySelector('#methods .panel p').textContent,'text/plain');
document.querySelectorAll('nav a').forEach(a=>a.onclick=()=>{document.querySelectorAll('nav a').forEach(n=>n.removeAttribute('aria-current'));a.setAttribute('aria-current','location')});
</script></main></div></body></html>'''


def write_report(root, output):
    bundle=collect(root);out=Path(output);out.mkdir(parents=True,exist_ok=True)
    (out/'results.json').write_text(json.dumps(bundle,indent=2)+'\n',encoding='utf-8')
    (out/'report.html').write_text(render_report(bundle),encoding='utf-8')
    (out/'methods.txt').write_text(methods_text(bundle)+'\n',encoding='utf-8')
