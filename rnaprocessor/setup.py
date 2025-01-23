from setuptools import setup, find_packages

setup(
    name="rnaprocessor",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        'biopython>=1.80',
        'numpy>=1.21.0',
        'pandas>=1.3.0',
        'matplotlib>=3.4.0',
        'seaborn>=0.11.0',
    ],
    entry_points={
        'console_scripts': [
            'rnaprocessor=rnaprocessor.__main__:main',
        ],
    },
    author="RNA Analysis Team",
    description="A tool for RNA sequence processing and BLAST analysis",
    python_requires=">=3.8",
) 