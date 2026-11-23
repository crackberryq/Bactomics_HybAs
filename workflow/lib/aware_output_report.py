"""Read-only, portable rendering of completed AWARE evidence."""
import csv
import html
import json
from pathlib import Path
from report_theme import evidence_page


def decisions(path, isolate):
    with Path(path).open() as stream:
        rows = list(csv.DictReader(stream, delimiter='\t'))
    if any(set(r) != {'field','value'} or r['value'] is None for r in rows):
        raise ValueError('Invalid AWARE filtering table')
    result = {r['field']: r['value'] for r in rows}
    if len(result) != len(rows) or result.get('isolate') != isolate:
        raise ValueError('Repeated filtering field or mismatched sample')
    return result


def render(document):
    esc = lambda value: html.escape(str(value if value not in (None, '') else 'Not recorded'), quote=True)
    def table(headers, rows):
        return '<div class="scroll"><table><thead><tr>' + ''.join('<th scope="col">'+esc(h)+'</th>' for h in headers) + '</tr></thead><tbody>' + ''.join('<tr>'+''.join('<td>'+esc(v)+'</td>' for v in row)+'</tr>' for row in rows) + '</tbody></table></div>'
    def percent(value):
        return 'Not available' if value is None else f'{100*value:.2f}%'
    def interval(value):
        return 'Not available' if not value else ' – '.join(percent(x) for x in value)
    actual = document.get('filtering', {})
    recommendation = document.get('parent_recommendation') or {}
    primary = next((g for g in document['groups'] if g['id']=='primary'), {})
    selected = actual.get('selected_target_taxid')
    selected_name = actual.get('selected_target_name')
    body = '<p class="warn">This report records the completed AWARE stage. It does not establish that assembly or downstream analysis completed. Check the saved run status in Bactologix.</p>'
    body += '<section id="target"><h2>Selected target and actual filtering</h2>'
    body += table(['Evidence','Value'], [
        ['AWARE+ mode',actual.get('aware_mode', 'legacy')],
        ['Selection mode',actual.get('kraken_decision_mode')],
        ['Selected target',f'{selected_name or "Not recorded"} · TaxID {selected or "Not recorded"}'],
        ['Selection source',actual.get('selected_target_source')],
        ['Automatic recommendation',f'{recommendation.get("name") or "Not recorded"} · TaxID {recommendation.get("taxid") or "Not recorded"}'],
        ['Assignment support',actual.get('decision_confidence')],
        ['Platform evidence', 'Single-platform evidence; cross-platform agreement unavailable' if actual.get('platform_agreement_reason') == 'single_platform_not_assessed' else actual.get('platform_agreement_reason')]])
    body += '<p>Assignment support describes database-dependent read assignments, not the probability of correct identification. All HybAs 10 input modes, polishing paths and AWARE+ features remain experimental and require validation.</p>'
    body += table(['Platform','Cleaning action','Reason','Retained units'], [[p,actual.get(p+'_kraken_decision'),actual.get(p+'_kraken_decision_reason'),f'{n:,} '+('pairs' if p=='illumina' else 'reads')] for p,n in primary.get('units',{}).items()])
    body += '<p>Actions above are the recorded per-platform cleaning decisions. Automatic recommendations and statistical evidence are shown separately.</p></section>'
    body += '<section id="groups"><h2>Retained groups</h2><p>'+esc(document['limitation'])+'</p><p>The primary output can overlap species groups; these are not independent replicates. Shared and unclassified groups are retained for inspection/export.</p>'
    body += table(['Group','Category','Retained units','Child route','Availability'], [[g['name'],g['category'],', '.join(p+': '+str(n)+(' pairs' if p=='illumina' else ' reads') for p,n in g['units'].items()),g.get('mode'),'Selectable for a new assembly draft' if g['assembly_eligible'] else g['reason']] for g in document['groups']])
    libraries=[]
    for group in document['groups']:
        for sample in (group.get('manifest') or {}).get('samples',[]):
            for lib in sample['libraries']:
                libraries.append([group['name'],lib['library_id'],lib['platform'],lib.get('chemistry',lib.get('library_type')),lib.get('orientation','—'),', '.join(g['run_id']+'/'+g['read_group_id'] for g in lib['read_groups'])])
    body += '<details><summary>Library and read-group provenance</summary>'+table(['Group','Library','Platform','Chemistry / type','Orientation','Run / read group'],libraries)+'</details></section>'
    body += '<section id="statistics"><h2>Statistical evidence</h2>'
    stats=document.get('statistics')
    if not stats:
        body += '<p>Statistical evidence was not requested for this run.</p>'
    else:
        body += '<p>'+esc(stats['interpretation'])+'</p><p>Statistical target: TaxID '+esc(stats['target_taxid'])+' · '+esc(stats['target_scope'])+'</p>'
        body += table(['Platform','Target units / total','Observed support','95% Wilson interval'],[[p,f'{v["target_units"]} / {v["units"]}',percent(v['fraction']),interval(v['wilson_95'])] for p,v in stats['platforms'].items()])
        stability=stats.get('stability')
        if stability:
            body += '<h3>Selection stability</h3>'+table(['Measure','Result'],[['Same automatic target',percent(stability['same_target_fraction'])],['Replicates',stability['replicates']],['Seed',stability['seed']],['Method',stability['method']]])
            body += table(['Selected TaxID','Replicates'],sorted(stability['selections'].items()))
        body += '<details><summary>Evidence by library and read group</summary>'+table(['Platform','Library','Run / read group','Target / total','Observed support','95% Wilson interval'],[[r['platform'],r['library'],r['run']+'/'+r['read_group'],f'{r["target_units"]} / {r["units"]}',percent(r['fraction']),interval(r['wilson_95'])] for r in stats['libraries']])+'</details>'
    body += '</section><section id="next"><h2>Next steps</h2><ol><li>Open AWARE outputs in Bactologix to select eligible groups and register derived samples.</li><li>Create separate run drafts for those samples. Their retained libraries determine the route; reclassification is optional.</li><li>For a parent stopped after AWARE, use Continue to assembly. The app checks saved evidence and required environments first.</li></ol><p>To move outputs to another workstation, copy this entire AWARE folder and import its manifest.json.</p><p><a href="manifest.json">Bundle manifest</a> · <a href="filtering.tsv">Recorded filtering decisions</a>'
    if stats:body+=' · <a href="statistics.json">Statistical evidence JSON</a>'
    body += '</p></section>'
    return evidence_page("AWARE target and filtering", document["isolate"] + " · AWARE stage complete · Experimental", body,
                         [("target","Target and filtering"),("groups","Retained groups"),("statistics","Statistical evidence"),("next","Next steps")])
