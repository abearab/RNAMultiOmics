from setuptools import setup, find_packages
from src import __version__, __author__, __email__
from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    python_requires='>=3.9',
    name='MultiOmics',
    description="A multi-omics analytic framework to study RNA dynamics",
    long_description=long_description,
    long_description_content_type='text/markdown',
    license="MIT License",

    version=__version__,
    author=__author__,
    author_email=__email__,
    maintainer=__author__,
    maintainer_email=__email__,

    url='https://github.com/abearab/MultiOmics',
    packages=find_packages(include=['src']),

    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ]
)