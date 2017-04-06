from collections import OrderedDict


class CompareSequence():
	'''
	args:
		template: List that contains the template's name and sequence
		alignments: List of List of alignments from the fasta file 
		domains = OrderedDict of domains
	'''
	def __init__(self, template, alignments, domains):
		self.template_name = template[0].replace('\n','')
		self.template_sequence = template[1].replace('\n','')
		self.alignments = alignments
		self.domains = domains
		self.gaps = 0
		self.gap_adjust = 0
		self.frontadjust = self.adjustPlacement()
		self.data = []
		self.TemplateVSAlignment()

	def TemplateVSAlignment(self):
		for alignment in self.alignments:
			self.gaps = 0
			self.data.append(self.DomainedSequences(alignment[0], alignment[1].replace('\n','')))
			print "\n"


	def adjustPlacement(self):
		'''
		checks to see where the first non gap is in a template
		'''
		return (i for i,v in enumerate(list(self.template_sequence)) if v!='-').next()
	
	def DomainedSequences(self, alignment, Oalignment):
		'''
		Loops through domains to get the correct substrings of the sequences
		based on the domain's logical specifications
		ie: domain = [1,10]		 <--- 10 Char in domain
		returns "T__ABCDL__EMST" <--- 10 Char in substring 
		'''
		combined_data = [alignment]
		for domain in self.domains:
			self.gap_adjust = 0
			current = self.template_sequence[self.frontadjust+self.gaps+self.domains[domain][0]-1:
										 	 self.frontadjust+self.gaps+self.domains[domain][1]]
			currentgaps = current.count('-')
			template, start, end = self.accountGaps(self.template_sequence, self.domains[domain], self.domains[domain], currentgaps)
			alignment = Oalignment[start:end]
			self.gaps += self.gap_adjust
			combined_data.extend(self.Analyze(self.domains[domain], template.upper(), alignment.upper()))
		return combined_data

	def Analyze(self, range, template, alignment):
		match = 0
		pcount = 0
		gcount = 0
		bottom_pcount = 0
		bottom_gcount = 0
		pmatch = 0
		gmatch = 0
		gap_matches = 0 
		for char in xrange(len(template)):
			if template[char] == alignment[char]:
				match += 1
				if template[char] == "-":
					gap_matches += 1
				if template[char] == "P" and alignment[char] == "P":
					pmatch += 1
				if template[char] == "G" and alignment[char] == "G":
					gmatch += 1
				if template[char] == "P":
					pcount += 1
				if template[char] == "G":
					gcount += 1
				if alignment[char] == "P":
					bottom_pcount += 1
				if alignment[char] == "G":
					bottom_gcount += 1
			else:
				if template[char] == "P" and alignment[char] == "P":
					pmatch += 1
				if template[char] == "G" and alignment[char] == "G":
					gmatch += 1
				if template[char] == "P":
					pcount += 1
				if template[char] == "G":
					gcount += 1
				if alignment[char] == "P":
					bottom_pcount += 1
				if alignment[char] == "G":
					bottom_gcount += 1
		percentage = round(float(match)/float(len(template))*100,2)
		data = [range, percentage, pmatch, gmatch, pcount, gcount,
          bottom_pcount, bottom_gcount, pcount - pmatch,
          gcount - gmatch, bottom_pcount - pmatch, bottom_gcount - gmatch]
		return data

	def accountGaps(self, Osequence, Odomain, domain, currentgaps):
		'''
		Recursive function that increases the accounted domain range (domain) until 
		the sequence String equals the length stated by the logical domain range (Odomain).
		Increases self.gaps to adjust to the total gaps in the corrected substring
		to keep the front adjustment accurate.

		args:
			Osequence: String - Sequence to account for. Either template or model.
			Odomain: List - Constant of original domain
			domain:  List - Domain that increases per iteration
			currentgaps: Int - current gaps present in the substring
		return:
			sequence: String - Correct substring adjusted with gaps accounted for
		modifies:
			self.gaps: Increased based on total gaps present in the corrected substring
		'''
		start = self.frontadjust + self.gaps + domain[0]-1
		end = self.frontadjust + self.gaps + domain[1]+currentgaps
		sequence = Osequence[start:end]
		gaps = sequence.count('-')

		if len(sequence.replace('-','')) == Odomain[1]- Odomain[0]+1:
			self.gap_adjust = gaps
			return sequence, start, end
		else:
			domain = [domain[0], domain[1]+((Odomain[1]- Odomain[0]+1)-len(sequence.replace('-','')))]
			return self.accountGaps(Osequence,Odomain, domain, gaps)