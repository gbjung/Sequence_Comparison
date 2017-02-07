from collections import OrderedDict
import csv
import sys
import re


def fileopen(file):
	'''
	opens file and parses the template and model sequences
	into appropriate dicts
	'''

	template = {}
	model = {}
	segment_start = []
	segment_end = []
	current_segment = 0
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

			elif "segments" in line:
				for line in f:
					if "start" in line:
						segment_start = map(int, line[8:].replace('\n','').split(','))
						current_segment = segment_start[0]
					if "end" in line:
						segment_end = map(int, line[6:].replace('\n','').split(','))
						break
					
			elif "structure" in line:
				for line in f:
					if line != "\n":
						for char in list(line.replace('\n', '')):
								'''
								print current_segment, segment_start
								if char == '/':
									current_segment = segment_start.pop(0)
								else:
									'''
								if char != '/':
									template[current_segment] = char
									current_segment += 1
					else:
						break
			elif "sequence" in line:
				current_segment = segment_start[0]
				for line in f:
					if line != "\n":
						#model.extend(line.replace('\n', ''))
						for char in list(line.replace('\n', '')):
							if char != '/':
								model[current_segment] = char
								current_segment += 1
					else:
						break

	return template, model, domains

def check_identity(template, model, domain, domain_name):

	correct = 0
	deletion = 0
	insertion = 0
	domain_index = range(domain[0],domain[1]+1)
	hold = ""
	print domain[0], domain[1]
	#template_data += "<span class = '{}'>".format(domain_name)
	#model_data += "<span class = '{}'>".format(domain_name)
	for sequence in domain_index:
		if template[sequence] == model[sequence]:
			hold += "[{}]".format(template[sequence])
			correct += 1
		elif template[sequence] == '-' and model[sequence] != '-':
			hold += template[sequence]
			deletion += 1
		elif template[sequence] != '-' and model[sequence] == '-':
			hold += template[sequence]
			insertion +=1
		else:
			hold += template[sequence]
	print hold
	return correct, deletion, insertion

def calculate(template, model, domain, domain_name):
	model = model
	template = template
	length_template = domain[1]-domain[0]
	length_model = domain[1]-domain[0]
	equivalent_residues, deletion, insertion = check_identity(template, model, domain, domain_name)
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
		for y in range(len(data[0])):
			wr.writerow([x[y] for x in data])

def generate_html(html):
	html = str(html)
	text_file = open("tableoutput.html", "w")
	text_file.write(html)
	text_file.close

def join(template, model):
	n = 10
	tem = "".join(template.values())
	tem = [tem[i:i+n] for i in range(0, len(tem), n)]

	mod = "".join(model.values())
	mod = [mod[i:i+n] for i in range(0, len(mod), n)]

	tem = [x for y in (tem[i:i+5] + ['<br>'] * (i < len(tem) - 2) for
     i in xrange(0, len(tem), 5)) for x in y]

	mod = [x for y in (mod[i:i+5] + ['<br>'] * (i < len(mod) - 2) for
     i in xrange(0, len(mod), 5)) for x in y]
	print tem, "\n", "\n", mod
	return mod, tem
if __name__ == "__main__":
	template_data = "<p>"
	model_data = ""
	file = sys.argv[1]
	template, model, domains = fileopen(file)
	data = results_dict(domains, template, model)
	export_to_csv(data)
	mod, tem = join(template, model)
	generate_html(tem)
