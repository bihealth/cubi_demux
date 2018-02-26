#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""Installation driver (and development utility entry point) for demuxtool
"""

from itertools import chain
import os
import sys

from setuptools.command.install import install
from setuptools import setup, find_packages
import pip
from pip.req import parse_requirements

import versioneer

__author__ = 'Manuel Holtgrewe <manuel.holtgrewe@bihealth.de>'

# Enforce python version >=3.4
if sys.version_info < (3, 4):
    print("At least Python 3.4 is required.\n", file=sys.stderr)
    sys.exit(1)

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

# Use DRMAA/non-DRMAA requirements

reqs = parse_requirements('requirements/base.txt', session=pip.download.PipSession())
requirements = [str(ir.req) for ir in reqs]

def console_scripts_entry_points(names, module):
    """Yield entries for the 'console_scripts' entry points"""
    for name in names:
        if module == 'wrappers':
            yield 'cubi-{name} = cubi_wrappers.{module}.{name}.__main__:main'.format(
                module=module, name=name)
        elif module == 'tools':
            yield 'cubi-{name} = cubi_wrappers.{module}.{name}:main'.format(
                module=module, name=name)


def bash_scripts(names):
    """Return generator with bash scripts of the given ``names``"""
    return (os.path.join('scripts', name) for name in names)


# Implementation of crazy baking-in of package repository path...


class TplToYamlInstall(install):
    """A custom install that expands *.yaml.tpl files to *.yaml files"""

    def initialize_options(self):
        """Set default values for options."""
        super().initialize_options()
        if not os.environ.get('BCL2FASTQ_CHANNEL'):
            assert False, 'BCL2FASTQ_CHANNEL not set'

    def run(self):
        res = super().run()
        print(vars(self))
        # Replace "__BCL2FASTQ_CHANNEL__" in *.yaml.tpl files
        if not getattr(self, 'manifest_files'):
            return res
        for package, paths in self.manifest_files.items():
            for path in paths:
                if path.endswith('.yaml.tpl'):
                    newpath = path[:-len('.tpl')]
                    print('=> path', path)
                    print('=> newpath', newpath)
                    with open(path, 'rt') as inputf:
                        inputs = inputf.read()
                    outputs = inputs.replace(
                        '__BCL2FASTQ_CHANNEL__', os.environ.get('BCL2FASTQ_CHANNEL'))
                    with open(newpath, 'wt') as outputf:
                        outputf.write(outputs)
        return res


cmdclass = versioneer.get_cmdclass().copy()
# cmdclass['install'] = TplToYamlInstall

setup(
    name='cubi_demux',
    version=versioneer.get_version(),
    cmdclass=cmdclass,
    description="CUBI Demuxtool",
    long_description=readme + '\n\n' + history,
    author="Manuel Holtgrewe",
    author_email='manuel.holtgrewe@bihealth.de',
    url='https://github.com/bihealth/demuxtool',
    packages=find_packages(),
    package_dir={
        'cubi_demux': 'cubi_demux',
    },
    entry_points={
        'console_scripts': (
            (
                'cubi-demux=cubi_demux.apps.cubi_demux:main',
            ),
        ),
    },
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords='cubi_pipeline',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
    ],
)

