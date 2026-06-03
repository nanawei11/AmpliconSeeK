import io
import os
from setuptools import find_packages, setup


here = os.path.abspath(os.path.dirname(__file__))

DESCRIPTION = (
    "AmpliconSeeK: a Python toolkit for detecting amplified genomic "
    "structures and candidate extrachromosomal DNA from sequencing data"
)

try:
    with io.open(os.path.join(here, "README.md"), encoding="utf-8") as f:
        long_description = "\n" + f.read()
except FileNotFoundError:
    long_description = DESCRIPTION

NAME = "ask-ecdna"
EMAIL = "nanawei11@163.com"
URL = "https://github.com/nanawei11/AmpliconSeeK/"
AUTHOR = "Nana Wei"
VERSION = "0.1.1"

setup(
    name=NAME,
    version=VERSION,
    author=AUTHOR,
    author_email=EMAIL,
    license="MIT",
    description=DESCRIPTION,
    url=URL,
    long_description_content_type="text/markdown",
    long_description=long_description,
    packages=find_packages(include=["ask", "ask.*", "data"]),
    package_data={"data": ["*"]},
    include_package_data=True,
    python_requires=">=3.8",
    install_requires=[
        "numpy",
        "pandas",
        "pysam",
        "matplotlib",
        "scipy",
        "statsmodels",
        "seaborn",
        "scikit-learn",
    ],
    entry_points={
        "console_scripts": [
            "ask=ask.cli:ask_main",
            "ask-search=ask.cli:ecdna_search_main",
        ],
    },
)
