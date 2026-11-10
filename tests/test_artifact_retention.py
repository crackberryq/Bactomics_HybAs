"""Storage contracts, including real Snakemake deletion and retry behavior."""
import json
import os
from pathlib import Path
import subprocess
import sys
import pytest
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'workflow/lib'))
from artifact_retention import archive_paths
from artifact_retention import _archive_job_paths

def test_move_dedup_and_mutation_isolation(tmp_path):
    native=tmp_path/'sample/work/tool';native.mkdir(parents=True)
    product=native/'native.sam';product.write_bytes(b'original')
    inode=product.stat().st_ino
    archive=archive_paths(tmp_path,[native,product],rule='tool')
    receipt=json.loads((archive/'manifest.json').read_text())
    assert receipt['state']=='complete' and len(receipt['entries'])==1
    saved=archive/'files/sample/work/tool/native.sam'
    assert saved.read_bytes()==b'original' and saved.stat().st_ino==inode
    native.mkdir();product.write_bytes(b'new')
    assert saved.read_bytes()==b'original'

def test_noop_and_unowned_paths(tmp_path):
    base=tmp_path/'base';base.mkdir();own=base/'old';own.write_text('keep')
    external=tmp_path/'external';external.write_text('outside')
    assert archive_paths(base,[base/'missing'],rule='none') is None
    with pytest.raises(ValueError,match='unowned'):
        archive_paths(base,[own,external],rule='bad')
    with pytest.raises(ValueError,match='input'):
        archive_paths(base,[own],rule='bad',inputs=[own])
    assert own.read_text()=='keep' and external.read_text()=='outside'
    assert not (base/'retained').exists()

def test_symlink_ancestor_is_refused(tmp_path):
    external=tmp_path/'external';external.mkdir();(external/'data').write_text('untouched')
    base=tmp_path/'base';base.mkdir();(base/'link').symlink_to(external,target_is_directory=True)
    with pytest.raises(ValueError,match='symlinked'):
        archive_paths(base,[base/'link/data'],rule='bad')
    assert (external/'data').read_text()=='untouched'

def test_failed_move_retains_source(tmp_path,monkeypatch):
    product=tmp_path/'result';product.write_text('old');original=Path.rename
    def fail(self,target):
        if self==product: raise OSError('storage failure')
        return original(self,target)
    monkeypatch.setattr(Path,'rename',fail)
    with pytest.raises(RuntimeError,match='no replacement'):
        archive_paths(tmp_path,[product],rule='tool')
    assert product.read_text()=='old'
    manifest=next((tmp_path/'retained').glob('*/manifest.json'))
    assert json.loads(manifest.read_text())['state']=='interrupted'

def invoke(tmp_path,*args,ok=True):
    snake=Path(sys.executable).with_name('snakemake')
    result=subprocess.run([str(snake),'-s',str(tmp_path/'Snakefile'),'--cores','1','--latency-wait','1',*args],cwd=tmp_path,text=True,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,timeout=60)
    assert (result.returncode==0)==ok,result.stdout
    return result

def make_workflow(tmp_path,body):
    (tmp_path/'Snakefile').write_text('import sys\nfrom pathlib import Path\n'+f'sys.path.insert(0, {str(ROOT / "workflow/lib")!r})\n'+'from artifact_retention import install\ninstall(workflow,str(Path.cwd()))\n'+body)

def test_real_engine_archive_before_rerun_and_noop(tmp_path):
    make_workflow(tmp_path,'''\nrule produce:
    output: directory("native")
    log: "tool.log"
    params: retention_paths=["sidecar"]
    shell: "mkdir -p native; printf new > native/result; printf new > sidecar; printf new > {log}"
''')
    (tmp_path/'native').mkdir();(tmp_path/'native/result').write_text('old result')
    (tmp_path/'native/internal.bin').write_text('native intermediate')
    (tmp_path/'sidecar').write_text('old sidecar');(tmp_path/'tool.log').write_text('old log')
    invoke(tmp_path,'--forceall','--dryrun')
    assert not (tmp_path/'retained').exists()
    invoke(tmp_path,'--forceall')
    archives=list((tmp_path/'retained').glob('*/files'));assert len(archives)==1
    assert (archives[0]/'native/internal.bin').read_text()=='native intermediate'
    assert (archives[0]/'tool.log').read_text()=='old log'
    assert (archives[0]/'sidecar').read_text()=='old sidecar'
    assert (tmp_path/'native/result').read_text()=='new'
    invoke(tmp_path)
    assert len(list((tmp_path/'retained').glob('*/files')))==1

def test_real_engine_failed_attempt_same_invocation_retry(tmp_path):
    make_workflow(tmp_path,'''\nrule produce:
    output: directory("native")
    log: "tool.log"
    shell: "mkdir -p native; if test ! -f attempted; then printf partial > native/partial; printf failed > {log}; touch attempted; exit 17; fi; printf success > native/result; printf success > {log}"
''')
    invoke(tmp_path,'--retries','1')
    archives=list((tmp_path/'retained').glob('*/files'));assert len(archives)==1
    assert (archives[0]/'native/partial').read_text()=='partial'
    assert (archives[0]/'tool.log').read_text()=='failed'
    assert not (tmp_path/'native/partial').exists()
    assert (tmp_path/'native/result').read_text()=='success'

def test_real_engine_keeps_failed_products_and_temporary_outputs(tmp_path):
    make_workflow(tmp_path,'''\nrule all:
    input: "final"
rule intermediate:
    output: temp("intermediate")
    shell: "printf intermediate > {output}"
rule finish:
    input: "intermediate"
    output: "final"
    shell: "printf partial > {output}; exit 17"
''')
    invoke(tmp_path,ok=False)
    assert (tmp_path/'intermediate').read_text()=='intermediate'
    assert (tmp_path/'final').read_text()=='partial'

def test_benchmark_is_archived_with_products(tmp_path):
    from artifact_retention import _job_paths
    from types import SimpleNamespace
    benchmark=tmp_path/'measurements.tsv';benchmark.write_text('old measurement')
    job=SimpleNamespace(output=[],log=[],benchmark=benchmark,params=SimpleNamespace(retention_paths=[]))
    archive=archive_paths(tmp_path,_job_paths(job),rule='measured')
    assert (archive/'files/measurements.tsv').read_text()=='old measurement'

def test_real_engine_python_run_rule_preserves_previous_product(tmp_path):
    make_workflow(tmp_path,'''\nrule produce:
    output: "result.json"
    run:
        Path(output[0]).write_text('new')
''')
    (tmp_path/'result.json').write_text('old')
    invoke(tmp_path,'--forceall')
    saved=list((tmp_path/'retained').glob('*/files/result.json'))
    assert len(saved)==1 and saved[0].read_text()=='old'
    assert (tmp_path/'result.json').read_text()=='new'


def test_registered_external_root_keeps_archive_outside_dataset(tmp_path):
    base = tmp_path / 'run'; base.mkdir()
    databases = tmp_path / 'databases'
    marker = databases / 'reference' / '.complete'
    marker.parent.mkdir(parents=True); marker.write_text('old marker')
    log = base / 'download.log'; log.write_text('old log')
    _archive_job_paths(base, [databases], [marker, log], rule='download')
    saved = list((databases / 'retained').glob('*/files/reference/.complete'))
    assert len(saved) == 1 and saved[0].read_text() == 'old marker'
    assert list((base / 'retained').glob('*/files/download.log'))[0].read_text() == 'old log'
    assert not (marker.parent / 'retained').exists()


def test_external_unapproved_path_refuses_entire_job_before_any_moves(tmp_path):
    base = tmp_path / 'run'; base.mkdir()
    local = base / 'product'; local.write_text('local bytes')
    external = tmp_path / 'database' / '.complete'
    external.parent.mkdir(); external.write_text('external marker')
    with pytest.raises(ValueError, match='unowned'):
        _archive_job_paths(base, [], [local, external], rule='download')
    assert local.read_text() == 'local bytes' and external.read_text() == 'external marker'
    assert not (base / 'retained').exists()


def test_most_specific_registered_root_owns_nested_product(tmp_path):
    databases = tmp_path / 'databases'
    marker = databases / 'reference' / '.complete'
    marker.parent.mkdir(parents=True); marker.write_text('old marker')
    _archive_job_paths(tmp_path, [databases], [marker], rule='download')
    assert list((databases / 'retained').glob('*/files/reference/.complete'))
    assert not (tmp_path / 'retained').exists()


def test_all_registered_groups_are_validated_before_first_move(tmp_path):
    base = tmp_path / 'run'; base.mkdir()
    databases = tmp_path / 'longer' / 'databases'; databases.mkdir(parents=True)
    marker = databases / '.complete'; marker.write_text('old marker')
    product = base / 'product'; product.write_text('shared bytes')
    alias = tmp_path / 'shared'; os.link(product, alias)
    with pytest.raises(ValueError, match='hardlinked'):
        _archive_job_paths(base, [databases], [marker, product], rule='download')
    assert marker.read_text() == 'old marker'
    assert product.read_text() == alias.read_text() == 'shared bytes'
    assert not (base / 'retained').exists() and not (databases / 'retained').exists()


@pytest.mark.parametrize('nested', [False, True])
def test_hardlinked_file_or_nested_product_is_refused_without_moves(tmp_path, nested):
    base = tmp_path / 'run'; base.mkdir()
    original = tmp_path / 'original'; original.write_text('shared bytes')
    path = base / 'native' / 'subdirectory' / 'output' if nested else base / 'output'
    path.parent.mkdir(parents=True, exist_ok=True)
    os.link(original, path)
    selection = base / 'native' if nested else path
    with pytest.raises(ValueError, match='hardlinked'):
        archive_paths(base, [selection], rule='tool')
    assert original.read_text() == path.read_text() == 'shared bytes'
    assert not (base / 'retained').exists()


def test_recursive_hardlink_check_does_not_follow_symlink_targets(tmp_path):
    base = tmp_path / 'run'; base.mkdir()
    external = tmp_path / 'external'; external.mkdir()
    source = external / 'product'; source.write_text('shared external bytes')
    os.link(source, external / 'alias')
    native = base / 'native'; native.mkdir()
    (native / 'directory-link').symlink_to(external, target_is_directory=True)
    (native / 'file-link').symlink_to(source)
    archive = archive_paths(base, [native], rule='tool')
    assert (archive / 'files/native/directory-link').is_symlink()
    assert (archive / 'files/native/file-link').is_symlink()
    assert source.read_text() == 'shared external bytes'


def test_registered_root_cannot_swallow_another_root_or_archive(tmp_path):
    databases = tmp_path / 'parent' / 'databases'; databases.mkdir(parents=True)
    marker = databases / '.complete'; marker.write_text('keep')
    with pytest.raises(ValueError, match='another owned root'):
        _archive_job_paths(tmp_path, [databases], [tmp_path / 'parent'], rule='bad')
    archive_root = tmp_path / 'retained'
    with pytest.raises(ValueError, match='retention archive'):
        _archive_job_paths(tmp_path, [archive_root], [marker], rule='bad')
    assert marker.read_text() == 'keep'


def test_real_engine_external_marker_requires_explicit_registered_root(tmp_path):
    base = tmp_path / 'run'; base.mkdir()
    databases = tmp_path / 'databases'
    marker = databases / 'reference' / '.complete'
    marker.parent.mkdir(parents=True); marker.write_text('old marker')
    make_workflow(base, '\nrule download:\n'
                  f'    output: {str(marker)!r}\n'
                  '    log: "tool.log"\n'
                  '    shell: "printf new > {output}; printf new-log > {log}"\n')
    snakefile = base / 'Snakefile'
    unregistered = snakefile.read_text()
    registered = unregistered.replace('install(workflow,str(Path.cwd()))',
                                      f'install(workflow,str(Path.cwd()),owned_roots=({str(databases)!r},))')
    snakefile.write_text(registered)
    invoke(base, '--forceall')
    saved = list((databases / 'retained').glob('*/files/reference/.complete'))
    assert len(saved) == 1 and saved[0].read_text() == 'old marker'
    assert marker.read_text() == 'new'
    snakefile.write_text(unregistered)
    result = invoke(base, '--forceall', ok=False)
    assert 'unowned output path' in result.stdout
    assert marker.read_text() == 'new'
    assert (base / 'tool.log').read_text() == 'new-log'
    assert len(list((databases / 'retained').iterdir())) == 1
