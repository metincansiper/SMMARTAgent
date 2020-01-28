import os
import urllib

PAX_JAR_URL = 'https://sourceforge.net/projects/biopax/files/paxtools/paxtools-5.0.1.jar/download'
PAX_JAR_PATH = 'jar/paxtools.jar'

def dl_paxtools_jar():
    urllib.request.urlretrieve (PAX_JAR_URL, PAX_JAR_PATH)

def lazy_dl_paxtools_jar():
    if not os.path.exists(PAX_JAR_PATH):
        dl_paxtools_jar()

lazy_dl_paxtools_jar()

os.environ['CLASSPATH'] = PAX_JAR_PATH
