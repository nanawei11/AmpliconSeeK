from setuptools import find_packages, setup


setup(
    name="ask-ecdna",
    version="0.1.0",
    description="ASK de novo ecDNA detection and targeted ASK-search",
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
