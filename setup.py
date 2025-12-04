
# Legacy setup.py for environments that still invoke setup.py directly.
from setuptools import setup, find_packages

setup(
    name="msaligner",
    version="1.0.1",
    description="Translation-based multiple sequence alignment toolkit",
    author="Ian Joens",
    license="MIT",
    package_dir={"": "src"},
    packages=find_packages(where="src"),
    include_package_data=True,
    python_requires=">=3.9",
    install_requires=[
        "biopython",
        "numpy",
        "matplotlib",
    ],
    entry_points={
        "console_scripts": [
            "msaligner=msaligner.cli:main",
        ]
    },
)
