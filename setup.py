from setuptools import setup
import pathlib

here = pathlib.Path(__file__).parent.resolve()
long_description = (here / "README.md").read_text(encoding="utf-8")

with open("requirements.txt") as requires_file:
    requirements = requires_file.read().split("\n")

setup(
    name             = 'namib',
    version          = '1.0.0',
    description      = 'Python package for handling and displaying multiple posterior distributions.',
    author           = 'Vasco Gennari',
    author_email     = 'vasco.gennari@gmail.com',
    url              = 'https://github.com/vascogennari/namib',
    long_description = long_description,
    packages         = ['namib'],
    python_requires  = '>=3',
    install_requires = requirements,
    entry_points     = {'console_scripts': ['namib = namib.namib:main']}
)