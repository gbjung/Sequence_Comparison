import os
import csv
from collections import OrderedDict
import sys
import re
from seleniumscript import TMHMM
from pprint import pprint

start = 0

def openDomains(file, fasta):
  '''
  args:
    file = Domains.txt
    fasta = GPCRdb_alignment.fasta
  returns:
    od = Domains from Domains.txt
    alignments = all the sequences from fasta file
    template = the template sequence from the fasta file from the Domain.txt specified,
               input for selenium script to get those domains
  '''
  od = ordered_dict() 
  alignments = []
  template = [] 
  sequence = "" #the template specified from the Domains.txt
  with open(file, 'rU') as f:
    for line in f:
      if "Sequence" in line:
        sequence = line[9:].strip()
      if "Domain" in line:
        domain_name = line.split(":")[0]
        domain_ranges = [x.replace("\n","").strip() for x in line.split(":")[1].split(',')]
        if len(domain_ranges) > 1:
          for drange in range(len(domain_ranges)):
            if "-" in domain_ranges[drange]:
              domain_ranges[drange] = [int(x) for x in domain_ranges[drange].split("-")]
              domain_ranges[drange] = [domain_ranges[drange][0], domain_ranges[drange][1]]
            else:
              domain_ranges[drange] = int(domain_ranges[drange])
        else:
          domain_ranges = [int(x) for x in domain_ranges[0].split("-")]
          domain_ranges = [domain_ranges[0], domain_ranges[1]]
        od[domain_name] = domain_ranges

    with open("GPCRdb_alignment.fasta", 'rU') as f:
      read = f.readlines()
      for line in range(len(read)):
        if sequence in read[line]:
          template = [read[line], read[line+1][::-1]]
          alignments.append([read[line][1:-1],read[line+1][1:-1]])
        elif ">" in read[line]:
          alignments.append([read[line][1:-1],read[line+1][1:-1]])
    return od, alignments, template


def ungapped(gapped_template, gapped_model, domain):
  template = gapped_template
  model = gapped_model 
  tem = []
  mod = []
  match = 0
  pcount = 0
  gcount = 0
  bottom_pcount = 0
  bottom_gcount = 0
  pmatch = 0
  gmatch = 0
  percentage = 0
  gaps_accounted = 0
  if len(domain)==2:
    tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch, gaps_accounted = main_algorithm(
      gapped_template, gapped_model, domain, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)

    percentage = round(float(match)/float(len(tem)-gaps_accounted)*100, 2)
  else:
    for num in domain:
      if isinstance(num,list):
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch, new_gaps_accounted = main_algorithm(
      gapped_template, gapped_model, num, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)
        gaps_accounted += new_gaps_accounted
      else:
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch = spot_alogrithm(num, gapped_template, gapped_model, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)
  
  print match
  print len(tem)
  print gaps_accounted

  percentage = round(float(match)/float(len(tem)-gaps_accounted)*100, 2)
  data = ["Match Percentage: {}".format(percentage),
          "Template: {}".format("".join(tem)), "Model: {}".format("".join(mod))]
          # "Template P Count: {}".format(pcount), "Template G Count: {}".format(gcount),
          # "Model P Count: {}".format(bottom_pcount), "Model G Count: {}".format(bottom_gcount),
          # "Template P Nonmatch: {}".format(pcount - pmatch), "Template G Nonmatch: {}".format(gcount - gmatch), 
          # "Model P Nonmatch: {}".format(bottom_pcount - pmatch), "Model G Nonmatch: {}".format(bottom_gcount - gmatch), 
          # ]

  return data

def gapped(template, model, domain):
  tem = []
  mod = []
  match = 0
  pcount = 0
  gcount = 0
  bottom_pcount = 0
  bottom_gcount = 0
  pmatch = 0
  gmatch = 0
  if len(domain)==2:
    tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch, gaps_accounted = main_algorithm(
      template, model, domain, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)

    percentage = round(float(match)/float(len(tem)-gaps_accounted)*100, 2)
  else:
    for num in domain:
      if isinstance(num,list):
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch, new_gaps_accounted = main_algorithm(
     template, model, num, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)
        gaps_accounted += new_gaps_accounted
      else:
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch = spot_alogrithm(num, template, model, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)
  
  percentage = round(float(match)/float(len(tem))*100, 2)
  data = ["Match Percentage: {}".format(percentage),
          "P Match: {}".format(pmatch), "G Match: {}".format(gmatch),
          "Template P Count: {}".format(pcount), "Template G Count: {}".format(gcount),
          "Model P Count: {}".format(bottom_pcount), "Model G Count: {}".format(bottom_gcount),
          "Template P Nonmatch: {}".format(pcount - pmatch), "Template G Nonmatch: {}".format(gcount - gmatch), 
          "Model P Nonmatch: {}".format(bottom_pcount - pmatch), "Model G Nonmatch: {}".format(bottom_gcount - gmatch), 
          ]

  return data

#checks to see if there are more gaps after the gaps have been accounted for
def checkgaps(template, domain, gaps, start):
  new_gaps = template[domain[0] - start + gaps:domain[1] - start + gaps].count("-")
  if new_gaps-gaps > 0:
    return checkgaps(template, [domain[0],domain[1]+new_gaps-gaps], new_gaps, start)
  else:
    return(new_gaps)

#checks to see how many previous chars were gaps and properly accounts for them
def gapinsurance(template, domain, gaps, start):
  print gaps
  new_gaps = template[0:domain[0] - start + gaps].count("-")
  if new_gaps > gaps:
    gaps = new_gaps
    return gapinsurance(template, domain, gaps, start)
  else:
    return gaps

#adjusts the custom domain to properly account for previous gaps
def adjust_custom(x, template, past_gaps):
  currentgaps = template[0:x-start+past_gaps].count("-")
  if currentgaps>past_gaps:
    past_gaps = past_gaps + currentgaps
    return adjust_custom(x,template, past_gaps)
  else:
    return currentgaps

#main algorithm for the alignment for when the domain is a normal range like A = [20,30]
def main_algorithm(gapped_template, gapped_model, domain, tem, mod, 
                   match, pcount, gcount, bottom_gcount, bottom_pcount, 
                   gmatch, pmatch):
  template = gapped_template
  model = gapped_model 
  gaps = template[0:domain[0] - start].count("-")
  moregaps = (gapinsurance(template, domain, gaps, start))
  gaps_accounted = (checkgaps(template, domain, moregaps, start))
  compare = range(domain[0]+moregaps, domain[1]+moregaps+gaps_accounted+1)
  for num in compare:
    if ("-" not in template[num - start]) and ("-" not in model[num - start]):
      if template[num - start] == model[num - start]:
        tem.append("["+template[num - start]+"]")
        mod.append("["+model[num - start]+"]")
        match+=1
        if template[num - start] == 'P':
          pmatch += 1
          pcount += 1
        if template[num - start] == 'G':
          gmatch += 1
          gcount += 1
        if model[num - start] == 'P':
          bottom_pcount += 1
        if model[num - start] == 'G':
          bottom_gcount += 1
      else:
        tem.append(template[num - start])
        mod.append(model[num - start])
        if template[num - start] == 'P':
          pcount += 1
        if template[num - start] == 'G':
          gcount += 1
        if model[num - start] == 'P':
          bottom_pcount += 1
        if model[num - start] == 'G':
          bottom_gcount += 1
    else:
        tem.append(template[num - start])
        mod.append(model[num - start])
  return tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch, gaps_accounted

def spot_alogrithm(num, gapped_template, gapped_model, tem, mod, 
                   match, pcount, gcount, bottom_gcount, bottom_pcount, 
                   gmatch, pmatch):
  template = gapped_template
  model = gapped_model 
  num = num+adjust_custom(num,template,template[0:num-start].count("-"))-start
  if template[num] == model[num]:
    tem.append("["+template[num]+"]")
    mod.append("["+model[num]+"]")
    match+=1
    if template[num] == 'P':
      pmatch += 1
      pcount += 1
    if template[num] == 'G':
      gmatch += 1
      gcount += 1
    if model[num] == 'P':
      bottom_pcount += 1
    if model[num] == 'G':
      bottom_gcount += 1
  else:
    tem.append(template[num])
    mod.append(model[num])
    if template[num - start] == 'P':
      pcount += 1
    if template[num - start] == 'G':
      gcount += 1
    if model[num - start] == 'P':
      bottom_pcount += 1
    if model[num - start] == 'G':
      bottom_gcount += 1
  return tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch

class ordered_dict(dict):
    def __init__(self, *args, **kwargs):
        dict.__init__(self, *args, **kwargs)
        self._order = self.keys()

    def __setitem__(self, key, value):
        dict.__setitem__(self, key, value)
        if key in self._order:
            self._order.remove(key)
        self._order.append(key)

    def __delitem__(self, key):
        dict.__delitem__(self, key)
        self._order.remove(key)

    def order(self):
        return self._order[:]

    def ordered_items(self):
        return [(key,self[key]) for key in self._order]


def loopthroughdomains(domains, template, model):
  '''
  args:
    domains: dict of domains
    template: the template sequence from the fasta file from the Domain.txt specified
    model: the model sequence to compare the template to
  returns:
    
  '''
  master = [model[0]]
  for item in domains:
    master.extend([item, "Range: {}".format(domains[item])]+ungapped(template[1].strip('\n'), model[1], domains[item]))
  return master

def exportCSV(data):
  with open('output.csv', 'wb') as output:
    writer = csv.writer(output)
    writer.writerows(data)

def main(textfile, fastafile):
  os.system('clear')
  custom_domains, alignments, template = openDomains(textfile, fastafile)
  #online_domains = TMHMM("".join(template)).domains
  total_data = []
  print template
  print alignments[0]
  total_data.append(loopthroughdomains(custom_domains, alignments[0], alignments[0]))
  print total_data
  # for alignment in alignments:
  #   total_data.append(loopthroughdomains(custom_domains, template, alignment))
  #print total_data
  exportCSV(total_data)

main(sys.argv[1], sys.argv[2])

