from setuptools import setup

setup(
    name='cfMEX_tools',
    version='0.3.0',
    description='A cell-free DNA multi-feature extraction tools',
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
