#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.md') as readme_file:
    readme = readme_file.read()

requirements = ['cclib>=1.7']

setup_requirements = [ ]

test_requirements = requirements.append(['pytest'])

setup(
    author="Alex M. Maldonado",
    author_email='aalexmmaldonado@gmail.com',
    python_requires='>=3',
    classifiers=[
        'Development Status :: 4 - Beta',
        'License :: OSI Approved :: MIT License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3'
    ],
    description="JSONs for computational chemistry",
    install_requires=requirements,
    extras_require={},
    license="MIT license",
    long_description=readme,
    include_package_data=True,
    keywords='qcjson',
    name='qcjson',
    packages=find_packages(include=['qcjson', 'qcjson.*']),
    setup_requires=setup_requirements,
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/keithgroup/qcjson',
    version='0.2.0',
    zip_safe=False,
)
