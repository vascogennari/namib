from setuptools import setup, find_packages
import pathlib, re

# here = pathlib.Path(__file__).parent.resolve()
# long_description = (here / "README.md").read_text(encoding="utf-8")

# with open("requirements.txt") as requires_file:
#     requirements = requires_file.read().split("\n")

def find_version(path, varname="__version__"):
    """Parse the version metadata variable in the given file.
    """
    with open(path, 'r') as fobj:
        version_file = fobj.read()
    version_match = re.search(
        r"^{0} = ['\"]([^'\"]*)['\"]".format(varname),
        version_file,
        re.M,
    )
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")

HERE = pathlib.Path(__file__).parent
with open(HERE / "pypi_description.rst", encoding='utf-8') as f:
    long_description = f.read()

with open(HERE / "requirements.txt") as requires_file:
    requirements = requires_file.read().split("\n")

setup(
    name             = 'namib',
    version          = find_version(HERE / "namib" / "__init__.py"),
    description      = 'Python package for handling and displaying multiple posterior distributions.',
    author           = 'Vasco Gennari',
    author_email     = 'vasco.gennari@gmail.com',
    url              = 'https://github.com/vascogennari/namib',
    long_description = long_description,
    packages         = find_packages(),
    python_requires  = '>=3',
    install_requires = requirements,
    entry_points     = {'console_scripts': ['namib = namib.namib:main']}
)