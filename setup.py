import subprocess
import os
path = subprocess.Popen('pwd', stdout=subprocess.PIPE).communicate()[0].strip()
os.system("cd {}; curl -O https://chromedriver.storage.googleapis.com/2.27/chromedriver_mac64.zip".format(path))
os.system("unzip -a chromedriver_mac64.zip; rm chromedriver_mac64.zip")
#os.system("mv chromedriver /usr/local/bin")