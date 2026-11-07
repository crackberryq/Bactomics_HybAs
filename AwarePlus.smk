"""Experimental preservation entry point, separate from the assembly workflow."""
from pathlib import Path
import sys

SOURCE = Path(workflow.basedir)
sys.path.insert(0, str(SOURCE / 'workflow/lib'))
from identification_review import read_json
from aware_plus import VERSION, SCHEMA

if not config.get('aware_plus_review') or not config.get('aware_plus_output'):
    raise ValueError('Supply aware_plus_review and a new aware_plus_output directory')
REVIEW = str(Path(config['aware_plus_review']).absolute())
OUTPUT = str(Path(config['aware_plus_output']).absolute())
DOCUMENT = read_json(REVIEW)
if not DOCUMENT.get('platforms', {}).get('ont'):
    raise ValueError('AWARE+ preservation requires ONT long reads for this sample')

rule all:
    input:
        OUTPUT

rule preserve_read_groups:
    input:
        review=REVIEW,
        evidence=[r['path'] for r in DOCUMENT['files']],
        source=[str(SOURCE / 'workflow/lib/aware_plus.py'),
                str(SOURCE / 'workflow/lib/identification_review.py'),
                str(SOURCE / 'workflow/scripts/aware_plus_preserve.py')]
    output:
        bundle=protected(directory(OUTPUT))
    params:
        version=VERSION, schema=SCHEMA
    conda:
        str(SOURCE / 'workflow/envs/controller.yml')
    threads: 1
    resources:
        mem_mb=1024
    log:
        OUTPUT + '.preservation.log'
    benchmark:
        OUTPUT + '.preservation.benchmark.tsv'
    shell:
        'python {input.source[2]:q} --review {input.review:q} --output {output.bundle:q} > {log:q} 2>&1'
