#!/usr/bin/env python3
"""Regression: stage a full install and repair stale Python files by reinstalling."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import shutil
import tempfile

p=argparse.ArgumentParser(description=__doc__)
p.add_argument('--build-dir',required=True,type=Path)
p.add_argument('--python',required=True)
p.add_argument('--cmake',default='cmake')
a=p.parse_args()
site=subprocess.check_output([a.python,'-c','import sysconfig; print(sysconfig.get_path("platlib"))'],text=True).strip()
with tempfile.TemporaryDirectory(prefix='mdcraft-install-test-') as temporary:
    env=dict(os.environ,DESTDIR=temporary)
    def install():
        subprocess.run([a.cmake,'--install',str(a.build_dir.resolve())],env=env,check=True,
                       stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    install()
    package=Path(temporary+site)/'mdcraft'
    init=package/'__init__.py';expected=init.read_bytes()
    init.write_text('STALE_INSTALL = True\n')
    extensions=list((package/'solver/potential').glob('_mdcraft_solver_potential*'))
    assert len(extensions)==1,extensions
    extension=extensions[0];extension.unlink()
    install()
    assert init.read_bytes()==expected,'__init__.py was not updated'
    assert extension.exists(),'Python extension was not reinstalled'
    stale = Path(temporary) / 'stale-libraries'
    stale.mkdir()
    library = next((package / '.libs').glob('libmdcraft-tools.so*'))
    copy = stale / library.name
    shutil.copy2(library, copy)
    with copy.open('ab') as target:
        target.write(b'old-build-test-marker')
    bad_env = dict(env, LD_LIBRARY_PATH=str(stale))
    result = subprocess.run([a.python, str(a.build_dir.resolve() / 'verify_python_install.py')],
                            env=bad_env, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    assert result.returncode and 'Remove conflicting LD_LIBRARY_PATH' in result.stdout, result.stdout
    # A subdirectory-only install must fail before installing runtime libraries.
    result=subprocess.run([a.cmake,'-P',str(a.build_dir.resolve()/'mdcraft/solver/cmake_install.cmake')],
                          env=env,stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True)
    assert result.returncode and 'TOP-LEVEL' in result.stdout,result.stdout
    print(json.dumps(dict(reinstall='passed',stale_init_repaired=True,missing_extension_repaired=True,
                          partial_install_rejected=True,stale_core_rejected=True,site_packages=site)))
