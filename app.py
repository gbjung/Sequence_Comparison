from collections import OrderedDict
import csv
import sys
import re

def fileopen(file):
	'''
	opens file and parses the template and model sequences
	into appropriate lists
	'''

	template = []
	model = []
	domains = OrderedDict()
	with open(file) as f:
		for line in f:
			if "Domains" in line:
				for line in f:
					if "start" in line:
						domain_name = line.split('_')[0]
						domain_start = [int(s) for s in line.split() if s.isdigit()][0]
						domain_end = [int(s) for s in line.split() if s.isdigit()][1]
						domains[domain_name] = domain_start, domain_end
					else:
						break
			elif "structure" in line:
				for line in f:
					if "*" in line:
						template.extend(line[:-3].replace('\n', ''))
						break
					else:
						template.extend(line.replace('\n', ''))
			elif "sequence" in line:
				for line in f:
					if "*" in line:
						model.extend(line[:-3].replace('\n', ''))
						break
					else:
						model.extend(line.replace('\n', ''))
	print template, model, domains
	return template, model, domains

def check_identity(template, model, domain):
	correct = 0
	deletion = 0
	insertion = 0
	tem = template[domain[0]:domain[1]]
	mod = model[domain[0]:domain[1]]
	for sequence in range(domain[1]-domain[0]):
		if template[sequence+domain[0]] == model[sequence+domain[0]]:
			tem[sequence] = '['+tem[sequence]+']'
			mod[sequence] = '['+mod[sequence]+']'
			correct += 1
		if template[sequence+domain[0]] == '-' and model[sequence+domain[0]] != '-':
			deletion += 1
		if tem[sequence] != '-' and mod[sequence] == '-':
			insertion +=1

	return correct, deletion, insertion

def calculate(template, model, domain, domain_name):
	model = model
	template = template
	length_template = domain[1]-domain[0]
	length_model = domain[1]-domain[0]
	equivalent_residues, deletion, insertion = check_identity(template, model, domain)
	try:
		percent_identity = float(equivalent_residues)/float(length_model)*100
	except ZeroDivisionError:
		percent_identity = 0
	d = OrderedDict()
	d["Domain"] = domain_name
	d["Length template"] = length_template
	d["Length model"] = length_model
	d["Deletions"] = deletion
	d["Insertions"] = insertion
	d["Equivalent Residues"] = equivalent_residues
	d["Percent Identity"] = percent_identity
	return d

def results_dict(domains, template, model):
	data = [['Domain', 'Length template', 'Length model',
			 'Deletions', 'Insertions', 'Equivalent Residues', 'Percent Identity']]

	for domain in domains:
		dataline = []
		domain_name = domain
		domain_values = [domains[domain][0], domains[domain][1]]
		results = calculate(template, model, domain_values, domain_name)
		for result in results:
			dataline.append(results[result])
		data.append(dataline)
	return data

def export_to_csv(data):
	with open('output.csv', 'wb') as myfile:
		wr = csv.writer(myfile)
		for y in range(len(data)):
			wr.writerow(data[y])
		

if __name__ == "__main__":
	file = sys.argv[1]
	template, model, domains = fileopen(file)
	data = results_dict(domains, template, model)
	export_to_csv(data)
