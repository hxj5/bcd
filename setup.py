# setup.py


import sys
from codecs import open
from os import path
from setuptools import setup, find_packages



here = path.abspath(path.dirname(__file__))


# load configures
exec(open("./bcd/app.py").read())


# Get the long description from the relevant file
with open(path.join(here, "README.rst"), encoding='utf-8') as f:
    long_description = f.read()


reqs = [
    "anndata", 
    "intervaltree", "matplotlib", "numpy", 
    "pandas", "scipy", "seaborn"
]


setup(
    name = APP,

    # Versions should comply with PEP440.  For a discussion on single-sourcing
    # the version across setup.py and the project code, see
    # https://packaging.python.org/en/latest/single_source_version.html
    version = VERSION,

    description = "bcd - Benchmarking of CNA Detection in Single-cell and Spatial Transcriptomics",
    long_description = long_description,
    long_description_content_type = "text/markdown",

    
    # The project's main homepage.
    url = "https://github.com/hxj5/bcd",

    # Author details
    author = 'Xianjie Huang, James Qiao',

    # Choose your license
    license='Apache-2.0',

    # What does your project relate to?
    keywords=['CNA', 'Benchmarking', "single cell", "spatial transcriptomics"],

    # You can just specify the packages manually here if your project is
    # simple. Or you can use find_packages().
    packages = find_packages(),

    #entry_points={
    #    'console_scripts': [
    #    ],
    #},

    # Check Python version. Only pip >= 9.0.1 supports.
    # Ref: https://stackoverflow.com/questions/42238484/prevent-package-from-being-installed-on-old-python-versions/42792413
    python_requires = ">=3.7",

    # List run-time dependencies here.  These will be installed by pip when
    # your project is installed. For an analysis of "install_requires" vs pip's
    # requirements files see:
    # https://packaging.python.org/en/latest/requirements.html
    
    install_requires = reqs,

    py_modules = ['bcd']

    # buid the distribution: python setup.py sdist
    # upload to pypi: twine upload dist/...
)
