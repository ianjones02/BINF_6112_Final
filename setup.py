from setuptools import setup, find_packages

setup(
    name="msaligner",
    version="1.0.0",
    description="Translation-Based Multiple Sequence Aligner",
    author="Ian Jones",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    install_requires=[
        "numpy>=1.21.0",
        "pandas>=1.3.0",
        "matplotlib>=3.5.0",
        "biopython>=1.79",
    ],
    entry_points={
        "console_scripts": [
            "msaligner=msaligner.cli:main",
        ],
    },
    python_requires=">=3.8",
)
