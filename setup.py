#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as fh:
	long_description = fh.read()

setup(
    name = "qed",
    version = "1.0.1",
    author = "Hans De Winter",
    author_email = "hans.dewinter@uantwerpen.be",
    description = ("Python implementation of the QED descriptor (Quantitative Estimation of Druglikeness)"),
    long_description = long_description,
	long_description_content_type="text/markdown",
    url = "https://github.com/UAMCAntwerpen/qed",
	download_url = "https://github.com/UAMCAntwerpen/qed/archive/1.0.1.tar.gz",
    packages = find_packages(include=['qed']),
	keywords = ['uamc', 'qed', 'rdkit', 'cheminformatics'],
    classifiers = [
        "Development Status :: 5 - Production/Stable",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
        "License :: Freeware",
		"Operating System :: POSIX :: Linux",
		"Operating System :: MacOS :: MacOS X",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Chemistry",
		"Topic :: Software Development :: Libraries :: Python Modules",
    ],
	python_requires = ">=3.6",
)