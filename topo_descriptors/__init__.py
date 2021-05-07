""" Initializations """

from pkg_resources import DistributionNotFound, get_distribution
from yaconfigobject import Config

try:
    __version__ = get_distribution(__name__).version
except DistributionNotFound:
    # package is not installed
    __version__ = ""

__author__ = """Daniele Nerini"""
__email__ = "daniele.nerini@meteoswiss.ch"

CFG = Config(name="topo_descriptors.conf")
