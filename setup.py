from setuptools import find_packages, setup

setup(
    name="refseqtools",
    packages=find_packages(),
    entry_points={
        "console_scripts": [
            "refseqtools = refseqtools.cli:cli"
        ]
    }
)
