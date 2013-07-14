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
      url = "https://github.com/Lawouach/WebSocket-for-Python",
      download_url = "http://www.defuze.org/oss/ws4py/",
      py_modules = ['emc'],
#      packages = ["ws4py", "ws4py.client", "ws4py.server"],
      platforms = ["any"],
      license = 'MIT',
      long_description = "WebSocket library for Python",
      classifiers=[
          'Development Status :: 4 - Beta',
          'Topic :: Software Development :: Libraries :: Python Modules'
          ],
     )
