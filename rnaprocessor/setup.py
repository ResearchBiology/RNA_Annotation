from setuptools import setup, find_packages

setup(
    name="new_rnaprocessor",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'biopython>=1.80',
        'tqdm>=4.65.0',
    ],
    python_requires=">=3.8",
)
