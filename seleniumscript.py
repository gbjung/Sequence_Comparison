from selenium import webdriver
from selenium.webdriver.support.ui import Select
from selenium.webdriver.common.keys import Keys
from selenium.webdriver.common.action_chains import ActionChains
from selenium.common.exceptions import NoSuchElementException
from selenium.webdriver.support import expected_conditions as EC
from selenium.webdriver.support.ui import WebDriverWait
import os
import subprocess

class TMHMM:
	def __init__(self, sequence):
		self.path = subprocess.Popen('pwd', stdout=subprocess.PIPE).communicate()[0].strip()
		self.chromedriver = self.path+"/chromedriver"
		self.sequence = sequence
		self.driver = webdriver.Chrome(self.chromedriver)
		self.base_url = "http://www.cbs.dtu.dk/services/TMHMM/"
		self.domains = self.cleanDomains(self.getDomains())


	def getDomains(self):
		d = self.driver
		stuff = self.sequence
		d.get(self.base_url)
		d.find_elements_by_xpath('//textarea[@name="SEQ"]')[0].send_keys(stuff)
		d.find_element_by_css_selector("input[type='radio'][value='-noplot']").click()
		d.find_element_by_css_selector("input[type='submit']").click()
		WebDriverWait(d, 10).until(EC.title_contains("result"))
		inputElement = d.find_element_by_css_selector("pre")
		domains = inputElement.text
		d.quit()
		return domains

	def cleanDomains(self, domains):
		clean_domains = {}
		dms = [x.encode('ascii','ignore') for x in domains.split("\n")]
		dms = [" ".join(x.split()) for x in dms]
		for domain in range(len(dms)):
			if "#" in dms[domain]:
				continue
			else:
				dms[domain] = dms[domain].split(" ")
				clean_domains["Domain_{}".format(domain-5)] = [int(dms[domain][-1]), int(dms[domain][-2])]
		return clean_domains
