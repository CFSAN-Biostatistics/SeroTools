#!/usr/bin/env python3


try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup


with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read().replace('.. :changelog:', '')

requirements = [
]

test_requirements = [
    "pytest",
]

setup(
    name='serotools',
    version='0.1.0',
    description="This package serves as a toolkit and repository for the White-Kauffmann-Le Minor scheme for Salmonella serotyping, which defines nomenclature and antigenic factors for each recognized serovar. The scheme is made available in multiple formats, along with methods for querying, comparing, and predicting serovar names and antigenic formulae.",
    long_description=readme + '\n\n' + history,
    author="Joseph D. Baugher, Ph.D.",
    author_email='joseph.baugher@fda.hhs.gov',
    url='https://github.com/CFSAN-Biostatistics/serotools',
    packages=['serotools'],
    package_dir={'serotools': 'serotools'},
    include_package_data=True,
    install_requires=requirements,
    license="BSD",
    zip_safe=False,
    keywords=['bioinformatics', 'Salmonella', 'serovar', 'serotype', 'serotyping', 'White-Kauffmann-Le Minor'],
    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    entry_points={'console_scripts': ['serotools = serotools.cli:main']},
    setup_requires=["pytest-runner"],
    tests_require=test_requirements
)
