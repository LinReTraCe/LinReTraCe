#! /usr/bin/env python

from setuptools import setup

from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(name='linretrace',
      version='1.5',
      description='Linear Response Transport Centre',
      long_description=long_description,
      long_description_content_type='text/markdown',
      author='Matthias Pickem',
      author_email='matthias.pickem@gmail.com',
      license='GPLv3',
      packages=['postproc','scattering','structure','structure.generators','structure.symmetries','structure.symmetries.onedim','structure.symmetries.twodim','structure.symmetries.threedim',
                'scripts','scripts.ltb','scripts.lwann','scripts.linspect'],
      install_requires=['numpy>=1.14','scipy>=1.10','h5py>=2.7','matplotlib>=2.2','ase>=3.17','spglib>=1.16'],
      scripts=['laverage','lconfig','ldft','linspect','lprint','ltb','ltb-run','lscat','lwann']
     )
