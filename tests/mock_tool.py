#!/usr/bin/env python3
"""Contract stubs for workflow orchestration tests. Never used in production.

These are deliberately not biological-tool simulations. Marker bases expose
which stage supplied the final output; SAM fixtures retain a secondary record.
"""
import gzip
import json
import os
from pathlib import Path
import shutil
import sys

tool = Path(sys.argv[0]).name
args = sys.argv[1:]


def arg(name, default=None):
    return args[args.index(name)+1] if name in args else default


def write(path, data):
    path=Path(path); path.parent.mkdir(parents=True, exist_ok=True); path.write_text(data)


def extend(path, base):
    lines=Path(path).read_text().splitlines()
    return '\n'.join(lines[:-1]+[lines[-1]+base])+'\n'


def records(path):
    opener=gzip.open if str(path).endswith('.gz') else open
    if Path(path).stat().st_size == 0: return []
    with opener(path,'rt') as handle: lines=handle.read().splitlines()
    return [lines[i+1] for i in range(0,len(lines),4)]


if os.environ.get('HYBAS_TEST_CALLS'):
    with open(os.environ['HYBAS_TEST_CALLS'],'a') as handle:
        handle.write(json.dumps({'tool':tool,'args':args})+'\n')
if tool == os.environ.get('HYBAS_TEST_FAIL_TOOL') and os.environ.get('HYBAS_TEST_FAIL_ISOLATE','') in ' '.join(args):
    print('Intentional workflow-test failure',file=sys.stderr)
    sys.exit(17)

if tool == 'conda':
    if '--version' in args: print('conda test-stub')
    elif 'list' in args:
        for name in ['unicycler','racon','medaka','polypolish','samtools','bwa','quast','seqkit','multiqc']:
            print(name+' 1.0 fixture')
    else: print('name: fixture\ndependencies: []')
elif tool == 'pigz':
    files=[];i=0
    while i<len(args):
        if args[i]=='-p': i+=2;continue
        if not args[i].startswith('-'):files.append(args[i])
        i+=1
    data=Path(files[-1]).read_bytes() if files else sys.stdin.buffer.read()
    encoded=gzip.compress(data,mtime=0)
    if '-f' in args and '-c' not in args:
        Path(files[-1]+'.gz').write_bytes(encoded);Path(files[-1]).unlink()
    else: sys.stdout.buffer.write(encoded)
elif tool == 'fastqc':
    write(Path(arg('-o'))/'fixture_fastqc.html','<html>fixture</html>')
elif tool == 'NanoPlot':
    write(Path(arg('-o'))/'NanoStats.txt','fixture')
elif tool == 'fastp':
    for src,dst in [('-i','-o'),('-I','-O')]:
        Path(arg(dst)).parent.mkdir(parents=True,exist_ok=True);shutil.copyfile(arg(src),arg(dst))
    write(arg('--html'),'<html>fixture</html>');write(arg('--json'),'{}')
elif tool == 'filtlong':
    with gzip.open(args[-1],'rb') as handle: sys.stdout.buffer.write(handle.read())
elif tool == 'kraken2':
    write(arg('--report'),'100\t2\t2\tG\t11\tExample\n')
    print('C\tfixture\t11\t4\t11:4')
elif tool == 'extract_kraken_reads.py':
    for src,dst in [('-s1','-o'),('-s2','-o2'),('-s','-o')]:
        if src in args:
            with gzip.open(arg(src),'rt') as handle: write(arg(dst),handle.read())
elif tool == 'unicycler':
    assert '--min_fasta_length' in args and '--spades_options' in args
    write(Path(arg('-o'))/'assembly.fasta','>fixture\n'+'A'*1200+'\n')
elif tool == 'minimap2': print('fixture\t4\t0\t4\t+\tfixture\t1200\t0\t4\t4\t4\t60')
elif tool == 'racon': sys.stdout.write(extend(args[-1],'C'))
elif tool == 'medaka_consensus': write(Path(arg('-o'))/'consensus.fasta',extend(arg('-d'),'G'))
elif tool == 'bwa':
    if args[0]=='index':
        for ext in ['amb','ann','bwt','pac','sa']:write(args[1]+'.'+ext,'fixture')
    else:
        fastqs=[x for x in args if x.endswith('.fq.gz')]
        print('@HD\tVN:1.6\tSO:unsorted\n@SQ\tSN:fixture\tLN:1200')
        for i,path in enumerate(fastqs):
            name='read'+str(i)
            if os.environ.get('HYBAS_TEST_NO_ALIGNMENTS'):
                print(f'{name}\t4\t*\t0\t0\t*\t*\t0\t0\tACGT\tIIII')
            else:
                print(f'{name}\t0\tfixture\t1\t60\t4M\t*\t0\t0\tACGT\tIIII\tNM:i:0')
                print(f'{name}\t256\tfixture\t11\t0\t4M\t*\t0\t0\t*\t*\tNM:i:0')
elif tool == 'samtools':
    if args[0]=='faidx':write(args[1]+'.fai','fixture\t1200\t0\t0\t0\n')
    elif args[0]=='index':write(args[1]+'.bai','fixture')
    elif args[0]=='sort':
        source=args[-1] if '-o' not in args or args.index('-o')<len(args)-2 else args[args.index('-o')-1]
        # Identify the source file rather than options or output.
        source=next(x for x in args[1:] if Path(x).is_file() and x!=arg('-o'))
        text=Path(source).read_text()
        if '-n' in args:text=text.replace('SO:unsorted','SO:queryname')
        write(arg('-o'),text)
    elif args[0]=='view':
        text=Path(args[-1]).read_text()
        lines=text.splitlines()
        if '-c' in args:
            data=[x for x in lines if not x.startswith('@')]
            if '-F' in args:data=[x for x in data if not int(x.split('\t')[1]) & int(arg('-F'))]
            print(len(data))
        elif '-H' in args:print('\n'.join(x for x in lines if x.startswith('@')))
        else:sys.stdout.write(text)
    else:print('samtools fixture')
elif tool == 'polypolish':
    assert args[0]=='polish' and len(args)==4, 'Polypolish must receive separate mate SAMs'
    for path in args[2:]:
        assert '\t256\t' in Path(path).read_text(), 'secondary alignments were lost'
        assert '\t*\t*\tNM:i:0' in Path(path).read_text(), 'SEQ=* records were lost'
    sys.stdout.write(extend(args[1],'T'))
elif tool == 'quast.py':
    out=Path(arg('-o'))
    length=sum(len(x) for x in Path(args[-1]).read_text().splitlines() if not x.startswith('>'))
    write(out/'report.tsv',f'Assembly\tfixture\n# contigs\t1\nTotal length\t{length}\nN50\t{length}\nGC (%)\t0\nLargest contig\t{length}\n')
    write(out/'report.txt','fixture\n')
elif tool == 'seqkit':
    path=args[-1]
    if path.endswith(('.fasta','.fa')):
        seqs=[''.join(x for x in Path(path).read_text().splitlines() if not x.startswith('>'))]
        fmt='FASTA'
    else:seqs=records(path);fmt='FASTQ'
    lengths=[len(x) for x in seqs]
    print('file\tformat\ttype\tnum_seqs\tsum_len\tmin_len\tavg_len\tmax_len')
    print(f'{path}\t{fmt}\tDNA\t{len(seqs)}\t{sum(lengths)}\t{min(lengths,default=0)}\t{sum(lengths)/max(1,len(lengths))}\t{max(lengths,default=0)}')
elif tool == 'multiqc':
    metadata=next(Path(x) for x in args if x.endswith('/metadata'))
    assert (metadata/'multiqc_custom_content.tsv').is_file(), 'MultiQC started before custom summary'
    assert '-c' in args
    write(Path(arg('-o'))/arg('-n'),'<html>fixture MultiQC</html>')
elif tool == 'prokka':
    for ext in ['gff','faa']:write(Path(arg('--outdir'))/(arg('--prefix')+'.'+ext),'fixture\n')
elif tool == 'busco':
    assert '--offline' in args, 'Unexpected BUSCO download during normal run'
    write(Path(arg('--out_path'))/arg('-o')/'short_summary.fixture.txt','C:99.0%[S:99.0%,D:0.0%],F:0.0%,M:1.0%,n:100\n')
else:
    raise SystemExit('Unhandled workflow-test executable: '+tool)
