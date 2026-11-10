"""Preserve reviewed native Python operations inside one owned output root.

Launchers verify native source before installing these module-local proxies.
This is not a filesystem monitor: native-child scratch and unwritten bytes are
outside its scope. Launchers must validate the persistent failure latch at exit.
Overwrite copies are independent, checksum-verified and fsynced before writing.
If finalization fails after a removal rename, the pending event names the retained
original; its former pathname can already be absent. Never delete pending events.
"""
from __future__ import annotations
import builtins
from datetime import datetime, timezone
import hashlib
import json
import mmap
import os
from pathlib import Path
import shutil
import stat
import subprocess
import uuid


class RetentionError(RuntimeError):
    """Preservation failed; never downgrade this to a missing native output."""


def _absolute(value):
    return Path(os.path.abspath(os.fsdecode(os.fspath(value))))


def _write_json(path, value):
    pending = path.with_name(path.name + '.pending-' + uuid.uuid4().hex)
    with builtins.open(pending, 'x', encoding='utf-8') as handle:
        json.dump(value, handle, sort_keys=True, indent=2)
        handle.write('\n')
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(pending, path)
    fd = os.open(path.parent, os.O_RDONLY | getattr(os, 'O_DIRECTORY', 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


def _fsync_directory(path):
    fd = os.open(path, os.O_RDONLY | getattr(os, 'O_DIRECTORY', 0))
    try:
        os.fsync(fd)
    finally:
        os.close(fd)


class _Proxy:
    def __init__(self, original, replacements):
        self._original, self._replacements = original, replacements

    def __getattr__(self, name):
        if name in self._replacements:
            return self._replacements[name]
        return getattr(self._original, name)


class RetentionSession:
    @classmethod
    def create(cls, root, *, tool, metadata=None):
        root = _absolute(root)
        if not root.is_dir() or root.resolve() != root:
            raise RetentionError('Native retention needs an existing output root without symlinked ancestors')
        parent = root / 'native-retained'
        if parent.is_symlink():
            raise RetentionError('Native retention directory must not be a symlink')
        parent.mkdir(exist_ok=True)
        path = parent / (datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S.%fZ') + '-' + uuid.uuid4().hex)
        path.mkdir(mode=0o700)
        (path / 'events').mkdir()
        with builtins.open(path / 'failure.bin', 'xb') as handle:
            handle.write(b'0')
            handle.flush()
            os.fsync(handle.fileno())
        _write_json(path / 'session.json', {'schema': 'hybas-native-retention-v1', 'root': str(root),
            'tool': str(tool), 'created_at': datetime.now(timezone.utc).isoformat(),
            'root_device': root.stat().st_dev, 'root_inode': root.stat().st_ino, 'metadata': metadata or {}})
        _fsync_directory(parent)
        _fsync_directory(root)
        return cls.attach(path)

    @classmethod
    def attach(cls, path=None):
        path = _absolute(path or os.environ['HYBAS_NATIVE_RETENTION_SESSION'])
        if path.resolve() != path:
            raise RetentionError('Native retention session must not contain symlinked ancestors')
        data = json.loads((path / 'session.json').read_text(encoding='utf-8'))
        root = _absolute(data['root'])
        if (data.get('schema') != 'hybas-native-retention-v1' or root.resolve() != root
                or path.parent != root / 'native-retained' or root.stat().st_dev != data['root_device']
                or root.stat().st_ino != data['root_inode']):
            raise RetentionError('Native retention session identity changed')
        session = cls()
        session.root, session.path, session.data = root, path, data
        session.env = {'HYBAS_NATIVE_RETENTION_SESSION': str(path)}
        fd = os.open(path / 'failure.bin', os.O_RDWR | getattr(os, 'O_NOFOLLOW', 0))
        try:
            if os.fstat(fd).st_size != 1:
                raise RetentionError('Invalid native retention failure latch')
            session._latch = mmap.mmap(fd, 1)
        finally:
            os.close(fd)
        return session

    def close(self):
        self._latch.close()

    def _fail(self, message, error=None):
        # Preallocated before work. Shared mmap remains visible to the parent
        # even if disk-full prevents writing an explanatory receipt.
        self._latch[0:1] = b'1'
        try:
            self._latch.flush()
        except OSError:
            pass
        raise RetentionError(message) from error

    def _path(self, value, *, writable=False):
        try:
            return self._checked_path(value, writable=writable)
        except RetentionError:
            raise
        except (OSError, ValueError, TypeError) as error:
            self._fail('Native output path could not be verified', error)

    def _checked_path(self, value, *, writable=False):
        path = _absolute(value)
        if (path == self.root or not path.is_relative_to(self.root)
                or path.is_relative_to(self.root / 'native-retained')):
            self._fail('Native retention refused an unowned or protected path: ' + str(path))
        if path.parent.resolve() != path.parent or self.root.resolve() != self.root:
            self._fail('Native retention refused symlinked ancestors: ' + str(path))
        if self.root.stat().st_ino != self.data['root_inode'] or self.root.stat().st_dev != self.data['root_device']:
            self._fail('Native output root identity changed')
        ancestor = path.parent
        while not ancestor.exists():
            ancestor = ancestor.parent
        if ancestor.stat().st_dev != self.data['root_device']:
            self._fail('Native retention refused a destination on another filesystem: ' + str(path))
        if writable and path.is_symlink():
            self._fail('Native retention refused writing through a symlink: ' + str(path))
        if writable and os.path.lexists(path):
            info = path.lstat()
            if (info.st_dev != self.data['root_device']
                    or stat.S_ISREG(info.st_mode) and info.st_nlink != 1
                    or not (stat.S_ISREG(info.st_mode) or stat.S_ISDIR(info.st_mode))):
                self._fail('Native retention refused an unowned writable product: ' + str(path))
        return path

    def _check_tree(self, path):
        try:
            self._check_tree_entries(path)
        except RetentionError:
            raise
        except OSError as error:
            self._fail('Native output tree could not be verified: ' + str(path), error)

    def _check_tree_entries(self, path):
        pending = [path]
        while pending:
            current = pending.pop()
            info = current.lstat()
            if info.st_dev != self.data['root_device']:
                self._fail('Native retention refused a product on another filesystem: ' + str(current))
            if stat.S_ISREG(info.st_mode):
                if info.st_nlink != 1:
                    self._fail('Native retention refused a hardlinked product: ' + str(current))
            elif stat.S_ISDIR(info.st_mode):
                with os.scandir(current) as entries:
                    pending.extend(Path(item.path) for item in entries)
            elif not stat.S_ISLNK(info.st_mode):
                self._fail('Native retention refused a special file: ' + str(current))

    def _event(self, operation, **details):
        destination = self.path / 'events' / uuid.uuid4().hex
        destination.mkdir()
        value = {'schema': 'hybas-native-retention-event-v1', 'state': 'pending', 'operation': operation,
            'created_at': datetime.now(timezone.utc).isoformat(), 'pid': os.getpid(), **details}
        _write_json(destination / 'manifest.json', value)
        _fsync_directory(destination.parent)
        return destination, value

    def _inventory(self, path, *, content=False, sync=False):
        try:
            return self._collect_inventory(path, content=content, sync=sync)
        except RetentionError:
            raise
        except OSError as error:
            self._fail('Native output inventory could not be verified: ' + str(path), error)

    def _collect_inventory(self, path, *, content=False, sync=False):
        """Never follow links. Detect mutation during inspection/copying."""
        entries, identities = {}, {}
        pending = [path]
        directories = []
        while pending:
            current = pending.pop()
            before = current.lstat()
            relative = current.relative_to(path).as_posix()
            identity = lambda info: (info.st_dev, info.st_ino, info.st_mode,
                info.st_nlink, info.st_size, info.st_mtime_ns, info.st_ctime_ns)
            identities[relative] = identity(before)
            if stat.S_ISREG(before.st_mode):
                if before.st_nlink != 1 or before.st_dev != self.data['root_device']:
                    raise OSError('Archive product identity is not independently owned: ' + str(current))
                entry = {'kind': 'file', 'bytes': before.st_size}
                if content or sync:
                    fd = os.open(current, os.O_RDONLY | getattr(os, 'O_NOFOLLOW', 0))
                    try:
                        if identity(os.fstat(fd)) != identity(before):
                            raise OSError('Archive product changed before inspection: ' + str(current))
                        if content:
                            digest = hashlib.sha256()
                            while True:
                                block = os.read(fd, 1024 * 1024)
                                if not block:
                                    break
                                digest.update(block)
                            entry['sha256'] = digest.hexdigest()
                        if sync:
                            os.fsync(fd)
                        if identity(os.fstat(fd)) != identity(before):
                            raise OSError('Archive product changed during inspection: ' + str(current))
                    finally:
                        os.close(fd)
                entries[relative] = entry
            elif stat.S_ISLNK(before.st_mode):
                entries[relative] = {'kind': 'symlink', 'target': os.readlink(current)}
            elif stat.S_ISDIR(before.st_mode):
                entries[relative] = {'kind': 'directory'}
                directories.append(current)
                with os.scandir(current) as children:
                    pending.extend(Path(child.path) for child in children)
            else:
                raise OSError('Archive contains a special file: ' + str(current))
        if sync:
            for directory in reversed(directories):
                _fsync_directory(directory)
            _fsync_directory(path.parent)
        return entries, identities

    def _preserve(self, path, *, operation, move):
        try:
            self._check_tree(path)
            destination, value = self._event(operation, original=str(path),
                original_relative=path.relative_to(self.root).as_posix(),
                retention='move-before-removal' if move else 'copy-before-overwrite')
            relative = 'files/' + path.relative_to(self.root).as_posix()
            saved = destination / relative
            saved.parent.mkdir(parents=True)
            if path.is_symlink():
                value.update(kind='symlink', link_target=os.readlink(path), target_contents_copied=False)
            else:
                value['kind'] = 'directory' if path.is_dir() else 'file'
            value['archive'] = relative
            _write_json(destination / 'manifest.json', value)
            before_copy = None if move else self._inventory(path)[1]
            if move:
                os.rename(path, saved)
            elif path.is_dir():
                shutil.copytree(path, saved, symlinks=True)
            elif path.is_symlink():
                saved.symlink_to(os.readlink(path))
            else:
                shutil.copy2(path, saved, follow_symlinks=False)
            inventory, _ = self._inventory(saved, content=not move, sync=True)
            if move:
                _fsync_directory(path.parent)
            else:
                source_inventory, after_copy = self._inventory(path, content=True)
                if source_inventory != inventory or before_copy != after_copy:
                    raise OSError('Native product changed or archive differs during copying')
            value['inventory'] = inventory
            ancestor = saved.parent
            while ancestor.is_relative_to(destination):
                _fsync_directory(ancestor)
                if ancestor == destination:
                    break
                ancestor = ancestor.parent
            value['state'] = 'complete'
            _write_json(destination / 'manifest.json', value)
            return saved
        except RetentionError:
            raise
        except Exception as error:
            self._fail('Native output preservation failed before ' + operation + ': ' + str(path), error)

    def remove(self, path, *, dir_fd=None):
        if dir_fd is not None:
            self._fail('Native retention does not accept directory-relative descriptors')
        path = self._path(path)
        if not os.path.lexists(path):
            return os.remove(path)
        if path.is_dir() and not path.is_symlink():
            raise IsADirectoryError(str(path))
        self._preserve(path, operation='remove', move=True)

    unlink = remove

    def rmtree(self, path, ignore_errors=False, onerror=None, *, onexc=None, dir_fd=None):
        if dir_fd is not None:
            self._fail('Native retention does not accept directory-relative descriptors')
        path = self._path(path)
        if not os.path.lexists(path) or path.is_symlink() or not path.is_dir():
            kwargs = {'ignore_errors': ignore_errors, 'onerror': onerror}
            if onexc is not None:
                kwargs['onexc'] = onexc
            return shutil.rmtree(path, **kwargs)
        self._preserve(path, operation='rmtree', move=True)

    def _before_write(self, path, operation):
        path = self._path(path, writable=True)
        if os.path.lexists(path) and not path.is_dir():
            self._preserve(path, operation=operation, move=False)
        return path

    def open(self, file, mode='r', *args, **kwargs):
        if isinstance(file, int):
            return builtins.open(file, mode, *args, **kwargs)
        if 'w' in mode or '+' in mode and 'a' not in mode:
            self._before_write(file, 'open-' + mode)
        elif any(flag in mode for flag in 'ax+'):
            self._path(file, writable=True)
        return builtins.open(file, mode, *args, **kwargs)

    def _copy(self, function, src, dst, *args, **kwargs):
        target = _absolute(dst)
        if function is not shutil.copyfile and target.is_dir():
            target = target / Path(src).name
        if os.path.exists(src) and os.path.exists(target) and os.path.samefile(src, target):
            return function(src, dst, *args, **kwargs)
        self._before_write(target, function.__name__)
        return function(src, dst, *args, **kwargs)

    def copyfile(self, src, dst, *args, **kwargs):
        return self._copy(shutil.copyfile, src, dst, *args, **kwargs)

    def copy(self, src, dst, *args, **kwargs):
        return self._copy(shutil.copy, src, dst, *args, **kwargs)

    def copy2(self, src, dst, *args, **kwargs):
        return self._copy(shutil.copy2, src, dst, *args, **kwargs)

    def copytree(self, src, dst, *args, **kwargs):
        target = self._path(dst, writable=True)
        # dirs_exist_ok is the fifth optional argument in shutil.copytree.
        merge = kwargs.get('dirs_exist_ok', args[4] if len(args) > 4 else False)
        if merge and os.path.lexists(target):
            inventory, _ = self._inventory(target)
            if any(entry['kind'] == 'symlink' for entry in inventory.values()):
                self._fail('Native retention refused merging into a tree containing symlinks: ' + str(target))
            self._preserve(target, operation='copytree', move=False)
        return shutil.copytree(src, dst, *args, **kwargs)

    def _rename(self, function, src, dst, *args, **kwargs):
        if kwargs.get('src_dir_fd') is not None or kwargs.get('dst_dir_fd') is not None:
            self._fail('Native retention does not accept directory-relative descriptors')
        source, target = self._path(src), self._path(dst, writable=True)
        if not os.path.lexists(source) or source == target:
            return function(src, dst, *args, **kwargs)
        self._check_tree(source)
        if os.path.lexists(target):
            if source.is_dir() != target.is_dir() or target.is_dir() and any(target.iterdir()):
                return function(src, dst, *args, **kwargs)
            self._preserve(target, operation=function.__name__, move=False)
        return function(src, dst, *args, **kwargs)

    def rename(self, src, dst, *args, **kwargs):
        return self._rename(os.rename, src, dst, *args, **kwargs)

    def replace(self, src, dst, *args, **kwargs):
        return self._rename(os.replace, src, dst, *args, **kwargs)

    def move(self, src, dst, *args, **kwargs):
        source = self._path(src)
        if os.path.lexists(source):
            self._check_tree(source)
        return self.move_preserving_destination(src, dst, *args, **kwargs)

    def move_preserving_destination(self, src, dst, *args, **kwargs):
        """Only for separately audited external-generated-log relocation sites."""
        source, target = _absolute(src), self._path(dst, writable=True)
        if target.is_dir():
            target = self._path(target / source.name, writable=True)
            if os.path.lexists(target):
                return shutil.move(src, dst, *args, **kwargs)
        if os.path.lexists(source) and os.path.lexists(target) and source != target:
            self._before_write(target, 'move')
        return shutil.move(src, dst, *args, **kwargs)

    def capture_process(self, command, result):
        """Keep unique captured streams; return object and command remain unchanged."""
        try:
            stdout = getattr(result, 'stdout', getattr(result, 'output', None))
            stderr = getattr(result, 'stderr', None)
            command = [os.fsdecode(os.fspath(part)) for part in command] if isinstance(command, (list, tuple)) else str(command)
            destination, value = self._event('native-process-streams', command=command,
                returncode=getattr(result, 'returncode', None), outcome=type(result).__name__,
                partial=isinstance(result, subprocess.TimeoutExpired), streams={})
            for name, data in (('stdout', stdout), ('stderr', stderr)):
                if data is None:
                    value['streams'][name] = {'captured': False}
                    continue
                if isinstance(data, str):
                    raw, encoding = data.encode('utf-8'), 'utf-8-from-decoded-text'
                elif isinstance(data, bytes):
                    raw, encoding = data, 'bytes'
                else:
                    raise TypeError('Native process stream is neither text nor bytes')
                with builtins.open(destination / (name + '.bin'), 'xb') as handle:
                    handle.write(raw)
                    handle.flush()
                    os.fsync(handle.fileno())
                value['streams'][name] = {'captured': True, 'path': name + '.bin', 'bytes': len(raw),
                    'sha256': hashlib.sha256(raw).hexdigest(), 'encoding': encoding}
            value['state'] = 'complete'
            _write_json(destination / 'manifest.json', value)
            return destination
        except RetentionError:
            raise
        except Exception as error:
            self._fail('Native process stream preservation failed', error)

    def os_proxy(self):
        return _Proxy(os, {'remove': self.remove, 'unlink': self.unlink, 'rename': self.rename, 'replace': self.replace})

    def shutil_proxy(self):
        return _Proxy(shutil, {'rmtree': self.rmtree, 'copyfile': self.copyfile,
            'copy': self.copy, 'copy2': self.copy2, 'copytree': self.copytree, 'move': self.move})

    def validate(self):
        if self._latch[:] != b'0':
            raise RetentionError('Native output preservation failed; inspect ' + str(self.path))
        for event in (self.path / 'events').iterdir():
            try:
                data = json.loads((event / 'manifest.json').read_text(encoding='utf-8'))
                if data.get('state') != 'complete':
                    raise ValueError('incomplete event')
            except (OSError, ValueError) as error:
                self._fail('Native output preservation has an incomplete event: ' + str(event), error)
        return True
