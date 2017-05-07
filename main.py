import os
import csv
from collections import OrderedDict
import sys
import re
from Comparison import CompareSequence

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
  od = OrderedDict() 
  alignments = []
  template = []
  online = False 
  sequence = "" #the template specified from the Domains.txt
  compare = "" #the specific comparison, optional
  with open(file, 'rU') as f:
    for line in f:
      if "Online" in line:
        online = True
      if "Sequence" in line:
        sequence = line[9:].strip()
      if "Extract" in line:
        compare = line[9:].strip()
      if "Domain" in line:
        domain_name = line.split(":")[0]
        domain_ranges = [x.replace("\n","").strip() for x in line.split(":")[1].split(',')]
        if len(domain_ranges) > 1:
          for drange in range(len(domain_ranges)):
            if "-" in domain_ranges[drange]:
              domain_ranges[drange] = [int(x) for x in domain_ranges[drange].split("-")]
              domain_ranges[drange] = [domain_ranges[drange][0], domain_ranges[drange][1]]
            else:
              domain_ranges[drange] = [int(domain_ranges[drange]), int(domain_ranges[drange])]
        else:
          domain_ranges = [int(x) for x in domain_ranges[0].split("-")]
          domain_ranges = [domain_ranges[0], domain_ranges[1]]
        od[domain_name] = domain_ranges
  with open("GPCRdb_alignment.fasta", 'rU') as f:
    read = f.readlines()
    for line in range(len(read)):
      if sequence in read[line]:
        template = [read[line][1:-1],read[line+1][1:-1]]
      if compare:
        if compare in read[line]:
          alignments = [[read[line][1:-1],read[line+1][1:-1]]]
      else:
        if ">" in read[line]:
          alignments.append([read[line][1:-1],read[line+1][1:-1]])

  return od, alignments, template, online

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
    master.extend(["Range: {}".format(domains[item])]+gapped(template[1].strip('\n'), model[1], domains[item]))
  return master

def exportCSV(filename, data):
  with open(filename, 'wb') as output:
    writer = csv.writer(output)
    writer.writerows(data)

def writeFile(temp,mod):
  with open('extracted.fasta','wb') as output:
    output.write(">"+temp[0]+'\n'+temp[1]+'\n')
    output.write(">"+mod[0]+'\n'+mod[1])
  output.close()

def main(textfile, fastafile):
  data = ["Range", "Match Percentage", "P Match","G Match", "Template P Count",
          "Template G Count", "Model P Count", "Model G Count",
          "Template P Nonmatch", "Template G Nonmatch", "Model P Nonmatch",
          "Model G Nonmatch"]
  shortdata = data[:2]

  custom_domains, alignments, template, online = openDomains(textfile, fastafile)
  if online:
    from seleniumscript import TMHMM
    online_domains = TMHMM("".join(template)).domains
  else:
    online_domains = {}
  total_domains = OrderedDict()
  for key in online_domains.keys() + custom_domains.keys():
    val_conct = ""
    try:
        val_conct = online_domains[key]
    except KeyError:
        val_conct = custom_domains[key]

    total_domains[key] = val_conct
  total_data = []
  total_shortdata = []
  domains_list = []
  domains_list.extend(["COMBINED"]*len(data))
  shortdomains_list = []
  shortdomains_list.extend(["COMBINED"]*len(shortdata))
  for domain in total_domains:
    domains_list.extend([domain]*len(data))
    shortdomains_list.extend([domain]*len(shortdata))
  compare = CompareSequence(template, alignments, total_domains)
  total_data.extend(compare.data)
  total_shortdata.extend(compare.shortdata)
  data = data*(len(total_domains)+1)
  shortdata = shortdata*(len(total_domains)+1)
  data.insert(0,"Sequence")
  shortdata.insert(0,"Sequence")
  total_data.insert(0, data)
  total_shortdata.insert(0, shortdata)
  domains_list.insert(0, "")
  shortdomains_list.insert(0,"")
  total_data.insert(0, domains_list)
  total_shortdata.insert(0, shortdomains_list)
  if compare.cleanedTemp and compare.cleanedMod:
    writeFile(compare.cleanedTemp,compare.cleanedMod)
  exportCSV('output.csv',total_data)
  exportCSV('shortended.csv',total_shortdata)

main(sys.argv[1], sys.argv[2])
