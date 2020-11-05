import os
from os import path, environ
import urllib
import jnius_config


file_dir = path.dirname(path.abspath(__file__))
root_dir = path.dirname(path.dirname(file_dir))

print(root_dir)

PAX_JAR_URL = 'https://sourceforge.net/projects/biopax/files/paxtools/paxtools-5.0.1.jar/download'
PAX_JAR_PATH = path.join(root_dir, 'jar/paxtools.jar')

def dl_paxtools_jar():
    urllib.request.urlretrieve (PAX_JAR_URL, PAX_JAR_PATH)

def lazy_dl_paxtools_jar():
    if not path.exists(PAX_JAR_PATH):
        dl_paxtools_jar()

lazy_dl_paxtools_jar()

cp_existing = os.environ.get('CLASSPATH')


if cp_existing is not None:
    os.environ['CLASSPATH'] =  cp_existing + ':' + PAX_JAR_PATH
else:
    os.environ['CLASSPATH'] = PAX_JAR_PATH
