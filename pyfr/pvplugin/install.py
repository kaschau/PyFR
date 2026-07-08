# Installer for the PyFR ParaView reader plugin
#
# Locates a ParaView installation, checks that its bundled h5py can read
# files written by the local h5py and, if it can not, installs a
# matching modern h5py wheel into a per-user site directory which the
# plugin falls back to at runtime.  Finally the plugin itself is copied
# to ~/.pyfr/paraview with the location of the pyfr package embedded.

import glob
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
import zipfile

# Wheel platform tags to request, by (sys.platform, machine)
_WHEEL_PLATFORMS = {
    ('darwin', 'arm64'): ['macosx_11_0_arm64', 'macosx_12_0_arm64',
                          'macosx_13_0_arm64', 'macosx_14_0_arm64'],
    ('darwin', 'x86_64'): ['macosx_10_9_x86_64', 'macosx_11_0_x86_64',
                           'macosx_12_0_x86_64'],
    ('linux', 'x86_64'): ['manylinux_2_17_x86_64', 'manylinux2014_x86_64',
                          'manylinux_2_28_x86_64', 'manylinux_2_34_x86_64'],
    ('linux', 'aarch64'): ['manylinux_2_17_aarch64', 'manylinux2014_aarch64',
                           'manylinux_2_28_aarch64'],
    ('win32', 'AMD64'): ['win_amd64'],
}

_PROBE_SCRIPT = '''
import json, platform, sys
info = {'version': list(sys.version_info[:3]),
        'machine': platform.machine(), 'sysplat': sys.platform}
try:
    import h5py
    info['h5py'] = h5py.__version__
    info['hdf5'] = h5py.version.hdf5_version
except ImportError:
    info['h5py'] = info['hdf5'] = None
print('PYFRJSON:' + json.dumps(info))
'''

_READ_SCRIPT = '''
import sys
if {site!r}:
    sys.path.insert(0, {site!r})
import h5py
with h5py.File({probe!r}, 'r') as f:
    f['eles'][...]
print('PYFRREAD:OK:' + h5py.version.hdf5_version)
'''


def _search_pvpython():
    home = os.path.expanduser('~')

    if sys.platform == 'darwin':
        pats = [f'{d}/ParaView*.app/Contents/bin/pvpython'
                for d in ['/Applications', f'{home}/Applications']]
    elif sys.platform == 'win32':
        pats = [os.path.join(pf, 'ParaView*', 'bin', 'pvpython.exe')
                for v in ['PROGRAMFILES', 'PROGRAMFILES(X86)']
                if (pf := os.environ.get(v))]
    else:
        # Common locations for extracted ParaView release tarballs
        pats = [f'{d}/{n}*/bin/pvpython'
                for d in ['/opt', '/usr/local', home, f'{home}/.local',
                          f'{home}/opt']
                for n in ['ParaView', 'paraview']]

    return [f for p in pats for f in glob.glob(p)]


def find_pvpython(path=None):
    if path:
        for cand in [path,
                     os.path.join(path, 'Contents', 'bin', 'pvpython'),
                     os.path.join(path, 'bin', 'pvpython'),
                     os.path.join(path, 'bin', 'pvpython.exe')]:
            if os.path.isfile(cand) and os.access(cand, os.X_OK):
                return cand

        raise ValueError(f'Unable to locate pvpython under {path}')

    if (pvpython := shutil.which('pvpython')):
        return pvpython

    # Prefer the most recent version
    if (found := _search_pvpython()):
        return max(found, key=lambda p: [int(d) for d in re.findall(r'\d+',
                                                                    p)])

    raise ValueError('Unable to locate ParaView; pass --paraview')


def _run_pv(pvpython, script):
    return subprocess.run([pvpython, '-c', script], capture_output=True,
                          text=True, timeout=300)


def probe_paraview(pvpython):
    res = _run_pv(pvpython, _PROBE_SCRIPT)
    for line in res.stdout.splitlines():
        if line.startswith('PYFRJSON:'):
            return json.loads(line[9:])

    raise RuntimeError(f'Unable to query {pvpython}:\n{res.stderr}')


def write_probe_file(path):
    import h5py
    import numpy as np

    # Representative of the most complex dtype in a PyFR file, written
    # with the same library version bounds as the native writers
    dt = np.dtype([('nodes', 'i8', (10,)), ('curved', '?'),
                   ('faces', [('cidx', 'i2'), ('off', 'i8')], (4,)),
                   ('colour', 'u1'), ('tags', 'u8')])
    with h5py.File(path, 'w', libver='latest') as f:
        f.create_dataset('eles', data=np.zeros(4, dtype=dt))


def pv_can_read(pvpython, probe, site=None):
    res = _run_pv(pvpython, _READ_SCRIPT.format(site=site or '', probe=probe))
    return any(l.startswith('PYFRREAD:OK') for l in res.stdout.splitlines())


def site_dir(info):
    major, minor = info['version'][:2]
    return os.path.expanduser(f'~/.pyfr/pv-site-{major}.{minor}')


def install_h5py(info, sitedir):
    major, minor = info['version'][:2]
    key = (info['sysplat'], info['machine'])

    try:
        platforms = _WHEEL_PLATFORMS[key]
    except KeyError:
        raise RuntimeError(f'No known h5py wheel platforms for {key}; '
                           f'install h5py manually and set PYFR_PV_SITE')

    with tempfile.TemporaryDirectory() as tmp:
        cmd = [sys.executable, '-m', 'pip', 'download', 'h5py',
               '--only-binary=:all:', '--no-deps', '--quiet',
               '--python-version', f'{major}.{minor}', '-d', tmp]
        for p in platforms:
            cmd += ['--platform', p]

        res = subprocess.run(cmd, capture_output=True, text=True)
        if res.returncode != 0:
            raise RuntimeError(f'pip download of h5py failed:\n{res.stderr}')

        wheel, = glob.glob(os.path.join(tmp, 'h5py-*.whl'))

        os.makedirs(sitedir, exist_ok=True)
        with zipfile.ZipFile(wheel) as z:
            z.extractall(sitedir)

        return os.path.basename(wheel)


def install_plugin(dstdir, dev=False):
    src = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                       'PyFR.py')

    import pyfr
    root = os.path.dirname(os.path.dirname(os.path.abspath(pyfr.__file__)))

    os.makedirs(dstdir, exist_ok=True)
    dst = os.path.join(dstdir, 'PyFR.py')

    if os.path.lexists(dst):
        os.remove(dst)

    if dev:
        # Symlink to the source tree; edits are picked up on the next
        # ParaView launch with no need to reinstall
        os.symlink(src, dst)
    else:
        text = open(src).read()
        text = text.replace('PYFR_ROOT = None', f'PYFR_ROOT = {root!r}', 1)
        with open(dst, 'w') as f:
            f.write(text)

    # Clean up artefacts of older installs
    for stale in ['pyfr_reader.py', 'pyfr-root']:
        if os.path.isfile(os.path.join(dstdir, stale)):
            os.remove(os.path.join(dstdir, stale))

    return dst


def _qsettings_file():
    if sys.platform == 'win32':
        base = os.environ.get('APPDATA') or os.path.expanduser('~')
        return os.path.join(base, 'ParaView', 'ParaView.ini')

    return os.path.expanduser('~/.config/ParaView/ParaView.ini')


def _paraview_exe(pvpython):
    bindir = os.path.dirname(os.path.abspath(pvpython))

    if sys.platform == 'darwin' and bindir.endswith('/Contents/bin'):
        return os.path.normpath(os.path.join(bindir, '..', 'MacOS',
                                             'paraview'))
    elif sys.platform == 'win32':
        return os.path.join(bindir, 'paraview.exe')
    else:
        return os.path.join(bindir, 'paraview')


# Our plugin entries inside a QSettings-escaped PluginsList XML value,
# along with any preceding entry separator
_PLUGIN_ENTRY_RE = re.compile(
    r'(?:\\n  )?<Plugin name=\\"(?:PyFR|pyfr_reader)\\"[^>]*?/>'
)


def register_autoload(pvpython, plugin):
    """Mark the plugin for auto-load in the ParaView GUI settings.

    The GUI persists its plugin list as an escaped XML string in the
    [PluginsList] section of its QSettings INI file, keyed by the path
    of the paraview executable; this performs the same edit the Plugin
    Manager makes when the "Auto Load" box is ticked.
    """
    ini = _qsettings_file()
    key = ('Local%3A'
           + _paraview_exe(pvpython).replace(':', '%3A').replace('/', '\\'))
    entry = (f'<Plugin name=\\"PyFR\\" filename=\\"{plugin}\\" '
             f'auto_load=\\"1\\" delayed_load=\\"0\\" version=\\"\\" '
             f'description=\\"\\" />')

    try:
        with open(ini) as f:
            lines = f.read().splitlines(keepends=True)
    except FileNotFoundError:
        lines = []

    kpfx = key + '='
    for i, l in enumerate(lines):
        if l.startswith(kpfx):
            val = _PLUGIN_ENTRY_RE.sub('', l[len(kpfx):].rstrip('\n'))
            if '\\n</Plugins>' not in val:
                break

            val = val.replace('\\n</Plugins>',
                              f'\\n  {entry}\\n</Plugins>', 1)
            lines[i] = kpfx + val + '\n'
            break
    else:
        block = (f'{kpfx}"<?xml version=\\"1.0\\" ?>\\n<Plugins>\\n'
                 f'  {entry}\\n</Plugins>\\n"\n')

        for i, l in enumerate(lines):
            if l.strip() == '[PluginsList]':
                lines.insert(i + 1, block)
                break
        else:
            if lines and not lines[-1].endswith('\n'):
                lines[-1] += '\n'
            lines += ['\n[PluginsList]\n', block]

    os.makedirs(os.path.dirname(ini), exist_ok=True)
    with open(ini, 'w') as f:
        f.writelines(lines)

    return ini


def process_install(args):
    pvpython = find_pvpython(args.paraview)
    info = probe_paraview(pvpython)
    major, minor = info['version'][:2]

    print(f'Found ParaView: {pvpython}')
    print(f'  Python {major}.{minor}, h5py {info["h5py"]}, '
          f'HDF5 {info["hdf5"]}')

    # Always provision the override h5py; PyFR files may originate from
    # systems with an HDF5 newer than the one ParaView bundles, and the
    # plugin only falls back to the override for files which need it
    with tempfile.TemporaryDirectory() as tmp:
        probe = os.path.join(tmp, 'probe.pyfrm')
        write_probe_file(probe)

        sitedir = site_dir(info)
        if (os.path.isdir(os.path.join(sitedir, 'h5py'))
                and pv_can_read(pvpython, probe, site=sitedir)):
            print(f'Override h5py already present: {sitedir}')
        else:
            wheel = install_h5py(info, sitedir)
            print(f'Override h5py installed: {wheel} -> {sitedir}')

            if not pv_can_read(pvpython, probe, site=sitedir):
                raise RuntimeError(f'h5py installed to {sitedir} but '
                                   f'verification failed')

        if pv_can_read(pvpython, probe):
            print('Note: bundled h5py also reads locally-written files; '
                  'the override will only be used if a file requires it')

    # Copy (or, in dev mode, symlink) the plugin
    dev = getattr(args, 'dev', False)
    dst = install_plugin(os.path.expanduser('~/.pyfr/paraview'), dev=dev)
    print(f'Plugin {"symlinked" if dev else "installed"}: {dst}')

    # Register it for auto-load in the GUI
    ini = register_autoload(pvpython, dst)
    print(f'Registered for auto-load in the ParaView GUI ({ini});\n'
          f'restart ParaView if it is currently running')

    print(f'''
For pvpython/pvbatch scripts:
  export PV_PLUGIN_PATH={os.path.dirname(dst)}''')


def process_status(args):
    try:
        pvpython = find_pvpython(args.paraview)
    except ValueError as e:
        print(e)
        return

    info = probe_paraview(pvpython)
    major, minor = info['version'][:2]

    print(f'ParaView: {pvpython}')
    print(f'  Python {major}.{minor}, h5py {info["h5py"]}, '
          f'HDF5 {info["hdf5"]}')

    with tempfile.TemporaryDirectory() as tmp:
        probe = os.path.join(tmp, 'probe.pyfrm')
        write_probe_file(probe)
        compat = pv_can_read(pvpython, probe)

        sitedir = site_dir(info)
        have_site = os.path.isdir(os.path.join(sitedir, 'h5py'))
        if have_site:
            ok = pv_can_read(pvpython, probe, site=sitedir)
            site_state = 'present, verified' if ok else 'present, BROKEN'
        else:
            site_state = 'absent'

    print(f'Bundled h5py reads locally-written files: '
          f'{"yes" if compat else "NO"}')
    print(f'Override h5py at {sitedir}: {site_state}')

    dst = os.path.expanduser('~/.pyfr/paraview/PyFR.py')
    print(f'Plugin at {dst}: '
          f'{"present" if os.path.isfile(dst) else "absent"}')

    try:
        with open(_qsettings_file()) as f:
            autoload = f'filename=\\"{dst}\\" auto_load=\\"1\\"' in f.read()
    except FileNotFoundError:
        autoload = False

    print(f'GUI auto-load registered: {"yes" if autoload else "no"}')
