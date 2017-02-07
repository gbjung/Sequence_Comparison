import os


start = 18
TM1 = [36,76]
TM2 = [77,86]
custom_domain [2, 10, 20, 54]
#Parameters: Analyze P's and G's, how many there are in total of individually for both sequences. How many of them are matching

#S = GHHHG #Get matching G's, number of G's, and num of unmatching G's for both
#A = GHHHH

#gaps dont' count, a match is only with residues aliging with residues
#ungapped: divide but ignore gaps in count 1/10
#gapped: 1/11
#get PG count for both
#fix ungapped
#custom domains 

YAY = (
"C-PKAGRHNY--IFVMIPTLYSIIFVVGIFGNSLVVIV---IYFY"+
"M-KLKTVASVFLLNLALADLCFLL-TLPLWAVYTA-MEYRWPFGNYLCKIA"+
"SASVSFNLYASVFLLTCLSIDRYLAIVHPMKS-----RLRRTMLVAKVTC"+
"IIIWLLAGLA-SLPAIIHRNVFFI/----------ITVCAF---HYE/TLPI"+
"GLGLTKNILGFLFPFLIILTSYTLIWKA-------LKK/-------------"+
"NDDIFKIIMAIVLFFFFSWIPHQIFTFLDVLIQLG------IIRDCRIAD"+
"IVDTAMPITICIAYFNNCLNPLFYG--FLGKKFKRYFLQLL--------")

B2 = (
"DVTQERDE-V-WVVGMGIVMSLIVLAIVFGNVLVITAIAK---F"+
"ERLQTVTNYFITSLACADLVMGLAVVPFGAAHILM--KMWTFGNFWCEFW"+
"TSIDVLCVTASIETLCVIAVDRYFAITSPFKYQSLLT-----KNKARVII"+
"LMVWIVSGLTSFLPIQM-HWYRAT-HQEAINCYANETCCDFFTN-------Q"+
"AYAIASSIVSFYVPLVIMVFVYSRVFQEAKRQLQK----IDKSEGR/KFCLK"+
"EHKALKTLGIIMGTFTLCWLPFFIVNIVHVI----QDNLIR---------"+
"--KEVYILLNWIGYVNSGFNPLIYCRS---PDFRIAFQELLCLRRSSLK")

def correct_ungapped(gapped_template, gapped_model, domain):
  template = gapped_template #gapped_template.replace("-","")
  model = gapped_model #gapped_model.replace("-","")
  tem = []
  mod = []
  match = 0
  pcount = 0
  gcount = 0
  pmatch = 0
  gmatch = 0
  gaps = template[0:domain[0] - start].count("-")
  gaps_accounted = (checkgaps(template, domain, gaps, start))
  compare = range(domain[0]+gaps, domain[1]+gaps_accounted+1)

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
      else:
        tem.append(template[num - start])
        mod.append(model[num - start])
        if template[num - start] == 'P':
          pcount += 1
        if template[num - start] == 'G':
          gcount += 1
    else:
        tem.append("("+template[num - start]+")")
        mod.append("("+model[num - start]+")")
  print("Correct Ungapped?: Ignores gaps, but does not remove them") 
  print("------------------")   
  print("4YAY : "+"".join(tem))
  print("BETA2: "+"".join(mod)+"\n")
  print("P Count: {}".format(pcount))
  print("G Count: {}".format(gcount))
  print("P Match: {}".format(pmatch))
  print("G Match: {}\n".format(gmatch))

  return round(float(match)/float(len(tem))*100, 2)

def gapped(template, model, domain):
  tem = []
  mod = []
  match = 0
  gaps = template[0:domain[0] - start].count("-")
  gaps_accounted = (checkgaps(template, domain, gaps, start))
  compare = range(domain[0]+gaps, domain[1]+gaps_accounted+1)
  for num in compare:
    if template[num - start] == model[num - start] and (template[num - start] != '-'):
      tem.append("["+template[num - start]+"]")
      mod.append("["+model[num - start]+"]")
      match+=1
    else:
      tem.append(template[num - start])
      mod.append(model[num - start])
  
  print("Gapped: Takes gaps into account, but not in domain length")
  print("------------------")      
  print("4YAY : "+"".join(tem))
  print("BETA2: "+"".join(mod)+"\n")

  return round(float(match)/float(len(tem)-gaps_accounted)*100, 2)
  # print("Gapped Match: {}".format(float(match)/float(len(tem))*100))
  # print("Ungapped Match: {}\n\n".format(float(match)/float(len([x for x in tem if x != "-"]))*100))

def checkgaps(template, domain, gaps, start):
  new_gaps = template[domain[0] - start:domain[1] - start + gaps].count("-")
  if new_gaps-gaps > 0:
    return checkgaps(template, [domain[0],domain[1]+new_gaps-gaps], new_gaps, start)
  else:
    return(new_gaps)

os.system('clear')
match_ungapped =ungapped(YAY, B2, TM1)
match_gapped = gapped(YAY, B2, TM1)
match_cungapped = correct_ungapped(YAY, B2, TM1)
print("Ungapped : {}%\nGapped : {}%\nCorrect Ungapped :{}%\n".format(
  match_ungapped,match_gapped,match_cungapped))
match_gapped = correct_ungapped(YAY, B2, TM2)
print match_gapped
