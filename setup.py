from setuptools import setup, find_packages
from codecs import open
from os import path

setup(
	name = 'loopDB',
	version = '0.2.2',
	description='A module for creating and storing DNA parts for Loop Assembly',
	url='https://github.com/HaseloffLab/LoopDB',
	download_url='https://github.com/HaseloffLab/LoopDB/archive/0.2.2.tar.gz',
	author = 'Mihails Delmans',
	author_email='md656@cam.ac.uk',
	license = 'GPL',
	classifiers=[
		# How mature is this project? Common values are
		#   3 - Alpha
		#   4 - Beta
		#   5 - Production/Stable
		'Development Status :: 4 - Beta',

		# Indicate who your project is intended for
		'Intended Audience :: Developers',
		'Intended Audience :: Science/Research',
		'Topic :: Scientific/Engineering :: Bio-Informatics',

		# Pick your license as you wish (should match "license" above)
		 'License :: OSI Approved :: MIT License',

		# Specify the Python versions you support here. In particular, ensure
		# that you indicate whether you support Python 2, Python 3 or both.
		'Programming Language :: Python :: 2.7',
	],
	keywords='bioinfomratics parts synthetic biology loop assembly',
	packages=find_packages(),
	install_requires=['sqlalchemy', 'partsdb'],
)

