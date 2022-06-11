#! /usr/bin/env python

from distutils.core import setup

setup(name='linretrace',
      version='1.1',
      description='Linear Response Transport Centre',
      author='Matthias Pickem',
      author_email='matthias.pickem@gmail.com',
      license='GPLv3',
      packages=['postproc','scattering','structure', 'structure.symmetries', 'structure.symmetries.onedim', 'structure.symmetries.twodim', 'structure.symmetries.threedim'],
      install_requires=['numpy','scipy','h5py','matplotlib','ase','spglib','boltztrap2'],
      scripts=['lconfig','ldft','lprint','ltb','lwann']
     )
