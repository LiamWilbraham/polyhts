from distutils.core import setup
from setuptools import find_packages

setup(
  name = 'polyhts',
  packages = ['polyhts'],
  version = '0.1',
  license='MIT',
  description = 'Polymer electronic property screening package.',
  author = 'Liam Wilbraham',
  author_email = 'liam.wilbrahaml@glasgow.ac.uk',
  url = 'https://github.com/LiamWilbraham/polyhts',
  download_url = 'https://github.com/LiamWilbraham/polyhts/archive/v_01.tar.gz',
  keywords = ['cheminformatics', 'chemistry', 'materials'],
  classifiers=[
    'License :: OSI Approved :: MIT License',
    'Programming Language :: Python :: 3',
  ],
)