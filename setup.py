from setuptools import setup, find_packages
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

with open("requirements.txt") as requires_file:
    requirements = requires_file.read().split("\n")

setup(
    name             = 'PYTOP',
    version          = '1.0.0',
    description      = 'Python package for handling and displaying multiple posterior distributions.',
    author           = 'Vasco Gennari',
    author_email     = 'vasco.gennari@ligo.org',
    long_description = long_description,
    packages         = ['PYTOP'],
    url              = 'https://github.com/vascogennari/PYTOP',
    python_requires  = '>=3',
    install_requires = requirements,
)