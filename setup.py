#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = ["Click>=7.0", "biopython==1.76", "pysam==0.15.3"]

setup_requirements = []

test_requirements = []

setup(
    author="Miles Smith",
    author_email="miles-smith@omrf.org",
    python_requires=">=3.5",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
    ],
    description="filter lines in a sam file that do not have a matching entry in a fasta file",
    entry_points={"console_scripts": ["filter_sam=filter_sam.cli:main",],},
    install_requires=requirements,
    license="BSD license",
    long_description=readme + "\n\n" + history,
    include_package_data=True,
    keywords="filter_sam",
    name="filter_sam",
    packages=find_packages(include=["filter_sam", "filter_sam.*"]),
    setup_requires=setup_requirements,
    test_suite="tests",
    tests_require=test_requirements,
    url="https://gitlab.com/milothepsychic/filter_sam",
    version="0.1.2",
    zip_safe=False,
)
