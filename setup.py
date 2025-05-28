from setuptools import setup

setup(
    name='cfMEX_tools',
    version='0.1',
    description='A sample Python package',
    author='Songphon Sutthitthasakul',
    author_email='songphon_sutthittha@cmu.ac.th',
    packages=['cfMEX_tools'],
    install_requires=[
        'numpy',
        'pandas',
        'biopython',
        'pysam'
    ],
)
