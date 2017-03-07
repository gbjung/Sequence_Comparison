import os
import csv
from collections import OrderedDict
import sys
import re
from seleniumscript import TMHMM

start = 18
TM1 = [25,35]
custom_domain = [21,22,25,27,30,31,32,33,44,56]
super_custom_domain = [21,22,[25,35],33,[44,56]]

stuff = (">5ht1b_human\n"+
"-------------------------------------------------------------------MEEPGAQCAPPPPAGSETWVPQANLSSAPSQNCSAKDYIYQDSIS-------LPWKVLLVMLLAL-ITLATTLSNAFVIATVYRT------RKLH------TPANYLIASLAVTDLLVSILVMPISTMYTVT-G----RW-TL-----GQVVCDFWLSSDITCCTASILHLCVIALDRYWAITDAVEYSAKR--T-PKRAAVMIALVWVFSISISLPPFFW------RQAKAEEEVSE-------------CVV------NTD------HILYTV-YSTVGAFYFPTLLLIALYGRIYVEARSRILKQTP--NRTGKRLTRAQLITDSPGSTSSVTSINSR---------------------------------------------------------------------------------------------------------------------------------------------------------VPDVPSESGSPVYVNQVKVRVSDALLEKK----KLMAARERKATKTLGIILGAFIV-CWLPFFIISLVMPI-----CKDA-CWF-------HLAIFDFFTWLGYLNSLINPIIYTMSNEDFKQAFHKLI-----------RFKCTS-----------------------------------------------------------------------------------------------------------------------")


YAY = (
"C-P--KAGRHNY--IFVMIPTLYSIIFVVGIFGNSLVVIV---IYFY"+  
"M-KLKTVASVFLLNLALADLCFLL-TLPLWAVYTA-MEYRWPFGNYLCKIA"+
"SASVSFNLYASVFLLTCLSIDRYLAIVHPMKS-----RLRRTMLVAKVTC"+
"IIIWLLAGLA-SLPAIIHRNVFFI/----------ITVCAF---HYE/TLPI"+
"GLGLTKNILGFLFPFLIILTSYTLIWKA-------LKK/-------------"+
"NDDIFKIIMAIVLFFFFSWIPHQIFTFLDVLIQLG------IIRDCRIAD"+
"IVDTAMPITICIAYFNNCLNPLFYG--FLGKKFKRYFLQLL--------")

B2 = (
"DVTQERDE-V-WV-VGMGIVMSLIVLAIVFGNVLVITAIAK---F"+
"ERLQTVTNYFITSLACADLVMGLAVVPFGAAHILM--KMWTFGNFWCEFW"+
"TSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLT-----KNKARVII"+
"LMVWIVSGLTSFLPIQM-HWYRAT-HQEAINCYANETCCDFFTN-------Q"+
"AYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQK----IDKSEGR/KFCLK"+
"EHKALKTLGIIMGTFTLCWLPFFIVNIVHVI----QDNLIR---------"+
"--KEVYILLNWIGYVNSGFNPLIYCRS---PDFRIAFQELLCLRRSSLK")


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
  template = {} 
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
        if "5ht1b_human" in read[line]:
          template[read[line]] = read[line+1][::-1]
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
    # print("Ungapped: Now properly accounts for gaps in the front of the domain and in the domain") 
  else:
    # print("Custom Domain: Alignment with a domain of custom points")
    # print("Custom Domain: {}".format(domain))
    for num in domain:
      if isinstance(num,list):
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch, new_gaps_accounted = main_algorithm(
      gapped_template, gapped_model, num, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)
        gaps_accounted += new_gaps_accounted
      else:
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch = spot_alogrithm(num, gapped_template, gapped_model, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)

    percentage = round(float(match)/float(len(tem)-gaps_accounted)*100, 2)

  # print("------------------")   
  # print("4YAY : "+"".join(tem))
  # print("BETA2: "+"".join(mod)+"\n")
  # print("First Sequence P Count: {}".format(pcount))
  # print("First Sequence G Count: {}".format(gcount))
  # print("Second Sequence P Count: {}".format(bottom_pcount))
  # print("Second Sequence G Count: {}".format(bottom_gcount))
  # print("P Match: {}".format(pmatch))
  # print("G Match: {}\n".format(gmatch))
  data = OrderedDict()
  data["type"] = "ungapped"
  data["tem"] = "".join(tem)
  data["mod"] = "".join(mod)
  data["pcount"] = pcount
  data["gcount"] = gcount
  data["nonmatchTopP"] = pcount - pmatch
  data["nonmatchTopG"] = gcount - gmatch
  data["bottom_pcount"] = bottom_pcount
  data["bottom_gcount"] = bottom_gcount
  data["nonmatchBotP"] = bottom_pcount - pmatch
  data["nonmatchBotG"] = bottom_pcount - gmatch
  data["pmatch"] = pmatch
  data["gmatch"] = gmatch
  data["percentage"] = percentage
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
    # print("Ungapped: Now properly accounts for gaps in the front of the domain and in the domain") 
  else:
    # print("Custom Domain: Alignment with a domain of custom points")
    # print("Custom Domain: {}".format(domain))
    for num in domain:
      if isinstance(num,list):
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch, new_gaps_accounted = main_algorithm(
     template, model, num, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)
        gaps_accounted += new_gaps_accounted
      else:
        tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch = spot_alogrithm(num, template, model, tem, mod, match, pcount, gcount, bottom_gcount, bottom_pcount, gmatch, pmatch)
  
  percentage = round(float(match)/float(len(tem))*100, 2)
  # print("Gapped: Aligned gaps aren't matches. Matches divided by domain length including gaps")
  # print("------------------")      
  # print("4YAY : "+"".join(tem))
  # print("BETA2: "+"".join(mod)+"\n")
  # print("First Sequence P Count: {}".format(pcount))
  # print("First Sequence G Count: {}".format(gcount))
  # print("Second Sequence P Count: {}".format(bottom_pcount))
  # print("Second Sequence G Count: {}".format(bottom_gcount))
  # print("P Match: {}".format(pmatch))
  # print("G Match: {}\n".format(gmatch))
  data = OrderedDict()
  data["type"] = "gapped"
  data["tem"] = "".join(tem)
  data["mod"] = "".join(mod)
  data["pcount"] = pcount
  data["gcount"] = gcount
  data["nonmatchTopP"] = pcount - pmatch
  data["nonmatchTopG"] = gcount - gmatch
  data["bottom_pcount"] = bottom_pcount
  data["bottom_gcount"] = bottom_gcount
  data["nonmatchBotP"] = bottom_pcount - pmatch
  data["nonmatchBotG"] = bottom_pcount - gmatch
  data["pmatch"] = pmatch
  data["gmatch"] = gmatch
  data["percentage"] = percentage

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
        # tem.append("("+template[num - start]+")")
        # mod.append("("+model[num - start]+")")
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

def export_to_csv(data):
  with open('output.csv', 'wb') as myfile:
    wr = csv.DictWriter(myfile, delimiter='\t', fieldnames=data)
    wr.writeheader()

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


os.system('clear')
match_ungapped = ungapped(YAY, B2, TM1)
match_gapped = gapped(YAY, B2, TM1)
custom_ungapped = ungapped(YAY,B2, custom_domain)
super_custom_ungapped = ungapped(YAY,B2, super_custom_domain)

a = TMHMM(stuff)

custom_domains, alignments, template = openDomains(sys.argv[1], sys.argv[2])
#export_to_csv(match_gapped)
print alignments
for domain in custom_domains:
  try:
    print domain[1]
    print ungapped(template, alignments[1][1], domain[1])
    print "\n"
  except IndexError:
    print("Out of range")