# -*- coding: utf-8 -*-
import sys
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

setup(name = "emc",
      version = '2.0-beta1',
      description = "Effective mass calculator",
      maintainer = "Alexandr Fonari",
      maintainer_email = "firstname.lastname@gatech.edu",
      url = "https://github.com/alexandr-fonari/emc",
      download_url = "https://github.com/alexandr-fonari/emc",
      py_modules = ['emc'],
      platforms = ["any"],
      license = 'MIT',
      long_description = "Effective mass calculator",
      classifiers=[
          'Development Status :: 4 - Beta',
          'Topic :: Scientific/Engineering :: Chemistry',
          'Topic :: Scientific/Engineering :: Physics'
          ],
     )
