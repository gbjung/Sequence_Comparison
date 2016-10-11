import sys
import os

class AutoVivification(dict):
    """Implementation of perl's autovivification feature."""
    def __getitem__(self, item):
        try:
            return dict.__getitem__(self, item)
        except KeyError:
            value = self[item] = type(self)()
            return value

def opendir():
	'''
	Goes through the directory and checks for template and model file
	'''
	for file in os.listdir(os.getcwd()):
		if "template" in file:
			try:
				template = fileopen(file)
			except IOError:
				print file + " doesn't exist?"
			except ValueError:
				print file + " isn't formated for this script"
		elif "model" in file:
			try:
				model = fileopen(file)
			except IOError:
				print file + " doesn't exist?"
			except ValueError:
				print file + " isn't formated for this script"

	return template, model


def fileopen(file):
	'''
	opens file function
	'''
	opened_file = list((open(file, "r").read().replace('\n', '')))
	return opened_file


def prompt_continue():
	'''
	prompts if the user wants to add more domains
	'''
	answer = raw_input("Add more domains? Y/N  : ").upper()
	if answer == 'Y' or answer == 'YES':
		return True
	elif answer == 'N' or answer == 'NO':
		return False
	else:
		print "Must type 'Y' or 'N' "
		return prompt_continue()

def check_if_start_valid(domain_name, current_segment, max_segment):
	'''
	checks if the domain start position is valid
	'''
	domain_name = domain_name
	current_segment = current_segment
	max_segment = max_segment

	domain_start = input("\n\nlast segment location: {} \nwhere does {} start?:  ".format(current_segment, domain_name))
	if domain_start <= current_segment:
		os.system('clear')
		print "The domain start must be greater than the current segment position"
		return check_if_start_valid(domain_name, current_segment, max_segment)
	if domain_start > max_segment:
		os.system('clear')
		print "The domain is out of range!"
		return check_if_start_valid(domain_name, current_segment, max_segment) 
	else:
		return domain_start

def check_if_end_valid(domain_name, current_segment, max_segment):
	'''
	checks if the domain end position is valid
	'''
	domain_name = domain_name
	current_segment = current_segment
	max_segment = max_segment

	domain_start = input("\n\nlast segment location: {} \nwhere does {} end?:  ".format(current_segment, domain_name))
	if domain_start <= current_segment:
		os.system('clear')
		print "The domain start must be greater than the current segment position"
		return check_if_start_valid(domain_name, current_segment, max_segment)
	if domain_start >= max_segment:
		os.system('clear')
		print "The domain is out of range!"
		return check_if_start_valid(domain_name, current_segment, max_segment) 
	else:
		return domain_start

def define_domains(template):
	'''
	creates a dictionary of domains and their start & end positions for compariosn
	'''
	temp_dict = AutoVivification()
	add_more = True
	current_segment = 0
	max_segment = len(template)
	os.system('clear')
	print "insert domains, insert done in domain to finish \n\n"
	while add_more:
		domain_name = raw_input("define the name of the domain:  ").strip()
		domain_start = check_if_start_valid(domain_name, current_segment, max_segment)
		current_segment = domain_start
		domain_end = check_if_end_valid(domain_name, current_segment, max_segment)
		current_segment = domain_end
		temp_dict[domain_name]["start"] = domain_start
		temp_dict[domain_name]["end"] = domain_end
		os.system('clear')
		print "{} out of {} total segements allocated in domain".format(current_segment,max_segment)
		add_more = prompt_continue()
	
	return temp_dict

def compare_domains(template, model, domains):
	'''
	makes comparisons of the model and template based on constraints from the domain
	'''
	for domain in domains:
		start = domains[domain]["start"]
		end = domains[domain]["end"]
		template_set = template[start:end]
		model_set = model[start:end]

		print "for domain {}\n".format(domain)
		print "{} match out of {} \n".format(len(set(template_set) & set(model_set)),len(template_set))
		print "{} percent match\n\n\n".format(100 * float(len(set(template_set) & set(model_set)))/float(len(template_set)))

if __name__ == "__main__":
	template = []	
	model = []
	domains = {}
	template, model = opendir()
	a = define_domains(template)
	domains.update(a)
	compare_domains(template, model, domains)
	
