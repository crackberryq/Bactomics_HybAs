"""Synthetic filesystem tests; no native analysis tools are executed."""
from pathlib import Path
import json
import multiprocessing
import os
import shutil
import stat
import subprocess
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'workflow' / 'lib'))
import native_output_retention as retention
from native_output_retention import RetentionError, RetentionSession


@pytest.fixture
def session(tmp_path):
    root = tmp_path / 'owned'
    root.mkdir()
    item = RetentionSession.create(root, tool='synthetic', metadata={'source': 'test'})
    yield item
    item.close()


def events(session):
    return [(directory, json.loads((directory / 'manifest.json').read_text()))
            for directory in sorted((session.path / 'events').iterdir())]


def saved(session):
    return [directory / value['archive'] for directory, value in events(session)
            if 'archive' in value]


def test_delete_repeated_generations_keeps_bytes_and_native_absence(session):
    product = session.root / 'result.txt'
    for data in (b'first', b'second'):
        product.write_bytes(data)
        inode = product.stat().st_ino
        session.remove(product)
        assert not product.exists()
        archive = next(path for path in saved(session) if path.read_bytes() == data)
        assert archive.stat().st_ino == inode
    assert {path.read_bytes() for path in saved(session)} == {b'first', b'second'}
    assert session.validate()
    assert json.loads((session.path / 'session.json').read_text())['metadata'] == {'source': 'test'}


def test_directory_retains_symlinks_without_copying_targets(session, tmp_path):
    outside = tmp_path / 'outside'
    outside.write_bytes(b'external')
    tree = session.root / 'stage'
    (tree / 'nested').mkdir(parents=True)
    (tree / 'nested' / 'bytes').write_bytes(b'generated')
    (tree / 'link').symlink_to(outside)
    session.rmtree(tree)
    assert not tree.exists()
    archive = saved(session)[0]
    assert (archive / 'nested' / 'bytes').read_bytes() == b'generated'
    assert (archive / 'link').is_symlink()
    assert os.readlink(archive / 'link') == str(outside)
    assert outside.read_bytes() == b'external'
    direct = session.root / 'direct-link'
    direct.symlink_to(outside)
    session.unlink(direct)
    assert not os.path.lexists(direct)
    assert any(value.get('target_contents_copied') is False for _, value in events(session))
    assert session.validate()


def test_overwrite_is_independent_and_append_has_no_snapshot(session):
    product = session.root / 'log'
    product.write_bytes(b'original')
    with session.open(product, 'wb') as handle:
        handle.write(b'new')
    archive = saved(session)[0]
    assert archive.read_bytes() == b'original'
    assert archive.stat().st_ino != product.stat().st_ino
    assert archive.stat().st_nlink == 1
    with session.open(product, 'ab') as handle:
        handle.write(b'-tail')
    assert len(saved(session)) == 1
    with pytest.raises(FileExistsError):
        session.open(product, 'xb')
    with session.open(product, 'r+b') as handle:
        handle.write(b'X')
    assert {item.read_bytes() for item in saved(session)} == {b'original', b'new-tail'}
    assert product.read_bytes() == b'Xew-tail'
    assert session.validate()


def test_native_missing_and_proxy_delegation(session, tmp_path):
    missing = session.root / 'absent'
    with pytest.raises(FileNotFoundError):
        session.os_proxy().unlink(missing)
    session.shutil_proxy().rmtree(missing, ignore_errors=True)
    with pytest.raises(FileNotFoundError):
        session.rmtree(missing)
    assert session.os_proxy().path is os.path
    assert session.shutil_proxy().which is shutil.which
    external = tmp_path / 'external'
    external.write_text('readable')
    with session.open(external) as handle:
        assert handle.read() == 'readable'
    assert not events(session)
    assert session.validate()


@pytest.mark.parametrize('operation', ['outside-remove', 'outside-open', 'outside-append', 'root-rmtree', 'protected'])
def test_unsafe_paths_latch_even_when_error_swallowed(session, tmp_path, operation):
    external = tmp_path / 'external'
    external.write_text('unmodified')
    with pytest.raises(RetentionError):
        if operation == 'outside-remove':
            session.remove(external)
        elif operation == 'outside-open':
            session.open(external, 'w')
        elif operation == 'outside-append':
            session.open(external, 'a')
        elif operation == 'root-rmtree':
            session.rmtree(session.root, ignore_errors=True)
        else:
            session.remove(session.path / 'failure.bin')
    assert external.read_text() == 'unmodified'
    with pytest.raises(RetentionError):
        session.validate()


def test_symlink_ancestors_and_write_through_symlinks_refused(session, tmp_path):
    external = tmp_path / 'external'
    external.mkdir()
    original = external / 'result'
    original.write_text('unmodified')
    link = session.root / 'alias'
    link.symlink_to(external, target_is_directory=True)
    with pytest.raises(RetentionError):
        session.remove(link / 'result')
    direct = session.root / 'direct'
    direct.symlink_to(original)
    with pytest.raises(RetentionError):
        session.open(direct, 'w')
    assert original.read_text() == 'unmodified'


@pytest.mark.parametrize('kind', ['hardlink', 'fifo'])
def test_unowned_entries_refused_before_directory_move(session, kind):
    tree = session.root / 'stage'
    tree.mkdir()
    (tree / 'good').write_text('untouched')
    if kind == 'hardlink':
        os.link(tree / 'good', tree / 'hard')
    else:
        os.mkfifo(tree / 'fifo')
    with pytest.raises(RetentionError):
        session.rmtree(tree, ignore_errors=True)
    assert (tree / 'good').read_text() == 'untouched'
    assert not saved(session)


@pytest.mark.parametrize('operation', ['remove', 'open'])
def test_io_failure_preserves_original_and_latches(session, monkeypatch, operation):
    original = session.root / 'result'
    original.write_bytes(b'valuable')
    def fail(*args, **kwargs):
        raise OSError('simulated full device')
    monkeypatch.setattr(retention.os if operation == 'remove' else retention.shutil,
                        'rename' if operation == 'remove' else 'copy2', fail)
    with pytest.raises(RetentionError):
        session.remove(original) if operation == 'remove' else session.open(original, 'w')
    assert original.read_bytes() == b'valuable'
    with pytest.raises(RetentionError):
        session.validate()


def test_payload_fsync_failure_does_not_truncate_original(session, monkeypatch):
    original = session.root / 'result'
    original.write_bytes(b'valuable')
    real_fsync = os.fsync
    def fail_archive_payload(fd):
        # Receipt files remain writable; fail only the independent archived copy.
        info = os.fstat(fd)
        if stat.S_ISREG(info.st_mode) and any(path.stat().st_ino == info.st_ino for path in saved(session) if path.exists()):
            raise OSError('simulated delayed payload write failure')
        return real_fsync(fd)
    monkeypatch.setattr(retention.os, 'fsync', fail_archive_payload)
    with pytest.raises(RetentionError):
        session.open(original, 'w')
    assert original.read_bytes() == b'valuable'
    with pytest.raises(RetentionError):
        session.validate()


def test_manifest_failure_is_visible_to_attached_parent(session, monkeypatch):
    original = session.root / 'result'
    original.write_bytes(b'valuable')
    attached = RetentionSession.attach(session.path)
    def fail(*args, **kwargs):
        raise OSError('simulated full device')
    monkeypatch.setattr(retention, '_write_json', fail)
    try:
        with pytest.raises(RetentionError):
            session.open(original, 'w')
        assert original.read_bytes() == b'valuable'
        with pytest.raises(RetentionError):
            attached.validate()
    finally:
        attached.close()


def test_copy_corruption_is_detected_before_truncation(session, monkeypatch):
    original = session.root / 'result'
    original.write_bytes(b'valuable')
    real_copy = shutil.copy2
    def corrupt(src, dst, **kwargs):
        result = real_copy(src, dst, **kwargs)
        Path(dst).write_bytes(b'corrupt!')
        return result
    monkeypatch.setattr(retention.shutil, 'copy2', corrupt)
    with pytest.raises(RetentionError):
        session.open(original, 'w')
    assert original.read_bytes() == b'valuable'
    with pytest.raises(RetentionError):
        session.validate()


def test_copy_replace_and_reviewed_external_move_keep_old_destination(session, tmp_path):
    target = session.root / 'output'
    target.write_bytes(b'old-output')
    source = tmp_path / 'input'
    source.write_bytes(b'input')
    session.shutil_proxy().copyfile(source, target)
    replacement = session.root / 'replacement'
    replacement.write_bytes(b'replacement')
    session.os_proxy().replace(replacement, target)
    external_log = tmp_path / 'generated.log'
    external_log.write_bytes(b'log')
    session.move_preserving_destination(external_log, target)
    assert target.read_bytes() == b'log'
    assert not external_log.exists()
    assert {path.read_bytes() for path in saved(session)} == {b'old-output', b'input', b'replacement'}
    assert session.validate()


def test_copytree_merge_preserves_old_tree(session, tmp_path):
    target = session.root / 'tree'
    target.mkdir()
    (target / 'result').write_bytes(b'old')
    source = tmp_path / 'new-tree'
    source.mkdir()
    (source / 'result').write_bytes(b'new')
    session.shutil_proxy().copytree(source, target, dirs_exist_ok=True)
    assert (target / 'result').read_bytes() == b'new'
    assert (saved(session)[0] / 'result').read_bytes() == b'old'
    assert session.validate()


def test_copytree_merge_never_writes_through_existing_symlink(session, tmp_path):
    target = session.root / 'tree'
    target.mkdir()
    external = tmp_path / 'external'
    external.write_bytes(b'untouched')
    (target / 'result').symlink_to(external)
    source = tmp_path / 'new-tree'
    source.mkdir()
    (source / 'result').write_bytes(b'new')
    with pytest.raises(RetentionError):
        session.copytree(source, target, dirs_exist_ok=True)
    assert external.read_bytes() == b'untouched'
    with pytest.raises(RetentionError):
        session.validate()


def test_append_refuses_hardlinked_file_without_modifying_either_name(session):
    product = session.root / 'product'
    product.write_bytes(b'old')
    os.link(product, session.root / 'alias')
    with pytest.raises(RetentionError):
        session.open(product, 'ab')
    assert product.read_bytes() == b'old'


def test_source_change_during_copy_prevents_destructive_open(session, monkeypatch):
    original = session.root / 'result'
    original.write_bytes(b'old')
    real_copy = shutil.copy2
    def change_source(src, dst, **kwargs):
        result = real_copy(src, dst, **kwargs)
        Path(src).write_bytes(b'changed')
        return result
    monkeypatch.setattr(retention.shutil, 'copy2', change_source)
    with pytest.raises(RetentionError):
        session.open(original, 'w')
    assert original.read_bytes() == b'changed'
    assert saved(session)[0].read_bytes() == b'old'
    with pytest.raises(RetentionError):
        session.validate()


def test_failure_after_move_leaves_recoverable_pending_archive(session, monkeypatch):
    original = session.root / 'result'
    original.write_bytes(b'valuable')
    real_write = retention._write_json
    def fail_completion(path, value):
        if value.get('state') == 'complete':
            raise OSError('simulated receipt commit failure after rename')
        return real_write(path, value)
    monkeypatch.setattr(retention, '_write_json', fail_completion)
    with pytest.raises(RetentionError):
        session.remove(original)
    assert not original.exists()
    assert saved(session)[0].read_bytes() == b'valuable'
    assert events(session)[0][1]['state'] == 'pending'
    with pytest.raises(RetentionError):
        session.validate()


@pytest.mark.parametrize('outcome', ['success', 'failure', 'timeout'])
def test_unique_process_capture_preserves_binary_streams(session, outcome):
    command = ['native', '--test']
    raw, errors = b'\x00\xffstdout', b'\x00\xfeerror'
    if outcome == 'success':
        result = subprocess.CompletedProcess(command, 0, raw, errors)
    elif outcome == 'failure':
        result = subprocess.CalledProcessError(3, command, output=raw, stderr=errors)
    else:
        result = subprocess.TimeoutExpired(command, 2, output=raw, stderr=errors)
    first = session.capture_process(command, result)
    second = session.capture_process(command, result)
    assert first != second
    assert result.stdout == raw and result.stderr == errors
    assert (first / 'stdout.bin').read_bytes() == raw
    assert (first / 'stderr.bin').read_bytes() == errors
    manifest = json.loads((first / 'manifest.json').read_text())
    assert manifest['partial'] is (outcome == 'timeout')
    assert manifest['command'] == command
    assert session.validate()


def _spawn_capture(path, index, fail):
    session = RetentionSession.attach(path)
    try:
        if fail:
            try:
                session.remove(session.root.parent / 'outside')
            except RuntimeError:
                pass  # Native code swallowing a retention failure cannot hide it.
        else:
            session.capture_process(['synthetic', str(index)],
                subprocess.CompletedProcess(['synthetic'], 0, str(index).encode(), b'error'))
    finally:
        session.close()


def test_spawn_workers_have_unique_streams_and_shared_failure_latch(session):
    context = multiprocessing.get_context('spawn')
    workers = [context.Process(target=_spawn_capture, args=(str(session.path), number, False)) for number in range(4)]
    for worker in workers:
        worker.start()
    for worker in workers:
        worker.join(15)
        assert worker.exitcode == 0
    assert len(events(session)) == 4
    assert {directory.joinpath('stdout.bin').read_bytes() for directory, _ in events(session)} == {b'0', b'1', b'2', b'3'}
    assert session.validate()
    worker = context.Process(target=_spawn_capture, args=(str(session.path), 0, True))
    worker.start()
    worker.join(15)
    assert worker.exitcode == 0
    with pytest.raises(RetentionError):
        session.validate()


def test_incomplete_event_cannot_pass_validation(session):
    session._event('interrupted')
    with pytest.raises(RetentionError):
        session.validate()
