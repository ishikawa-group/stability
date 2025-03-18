from setuptools import setup, find_packages

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="stability",
    version="0.0.1",
    author="Ishikawa-lab@ScienceTokyo",
    author_email="ishikawa.a.ai@m.titech.ac.jp",
    description="Material stability analysis tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
    python_requires=">=3.8",
    install_requires=[
        "ase>=3.22.1",
        "numpy>=1.26.2",
        "scipy>=1.11.3",
        "matplotlib>=3.8.1",
        "pymatgen>=2023.11.12",
        "mp_api>=0.45.3",
    ],
)
