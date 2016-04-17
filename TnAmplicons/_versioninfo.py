# _versioninfo.py
#
# gets the version number from the package info
# checks it agains the github version
import sys
from pkg_resources import get_distribution, parse_version

try:
    _dist = get_distribution('TnAmplicons')
    version_num = _dist.version
except:
    version_num = 'Please install this project with setup.py'


version_master = "https://raw.githubusercontent.com/msettles/TnAmplicons/master/VERSION"
repo_master = "https://github.com/msettles/TnAmplicons"
version_develop = "https://raw.githubusercontent.com/msettles/TnAmplicons/develop/VERSION"
repo_develop = "https://github.com/msettles/TnAmplicons/tree/develop"

try:
    import urllib2
    github_version_num = urllib2.urlopen(version_master).readline().strip()
    if parse_version(github_version_num) > _dist.parsed_version:
        sys.stderr.write("A newer version (%s) of TnAmplicons is available at %s\n" % (github_version_num, repo_master))
    elif parse_version(github_version_num) < _dist.parsed_version:
        github_version_num = urllib2.urlopen(version_develop).readline().strip()
        if parse_version(github_version_num) > _dist.parsed_version:
            sys.stderr.write("A newer version (%s) of TnAmplicons is available at %s\n" % (github_version_num, repo_develop))
except:
    sys.stderr.write("Error retrieving github version_number\n")

__version__ = version_num
