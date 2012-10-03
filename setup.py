# -*- coding: utf-8 -*-

from setuptools import setup, find_packages

setup(name='transitlib',
      version='0.1',
      description="Star planet transit calculation package",
      author="Anatoli Vladev",
      author_email="avladev@gmail.com",
      url="http://github.com/acshu/transit-lib",
      packages=find_packages(),
      scripts=['transit','transit.bat'],
      platforms="All",
      license="GPL",
      
      classifiers=[
      'Development Status :: 3 - Alpha',
      'Environment :: Console',
      'Intended Audience :: Education',
      'Intended Audience :: Science/Research',
      'License :: OSI Approved :: GNU General Public License (GPL)',
      'Natural Language :: English',
      'Operating System :: Unix',
      'Operating System :: Microsoft',
      'Operating System :: Other OS',
      'Programming Language :: Python',
      'Programming Language :: Python :: 2.5',
      'Programming Language :: Python :: 2.6',
      'Programming Language :: Python :: 2.7',
      'Topic :: Scientific/Engineering :: Astronomy',
      'Topic :: Utilities'
      ]
    )