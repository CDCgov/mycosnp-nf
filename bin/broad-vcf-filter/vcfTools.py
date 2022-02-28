# python vcf tools

from __future__ import division

import re, sys

from scipy.stats import binom_test

class VcfRecord:
	"""
	
	"""
	def __init__(self, vcf_line):
		self.vcf_line = vcf_line.rstrip()
		
		fields = self.vcf_line.split('\t')
		self.chrom = fields[0]
		self.pos = fields[1]
		self.id = fields[2]
		self.ref = fields[3]
		self.alt = fields[4]
		self.qual = fields[5]
		self.filter = fields[6]
		self.info = fields[7]
		self.format = fields[8]
		
		self.genotypes = []
		self.vcf_annot = False
		for i in range(9,len(fields)):
			if fields[i]:
				gt_exp = re.compile(r"^[\d\.]")
				if re.search(gt_exp,fields[i]):
					self.genotypes.append(fields[i])
				else:
					self.vcf_annot = fields[i]
					
	def is_passing(self,caller):
		if self.filter == 'PASS':
			return True
		elif caller == 'GATK' and self.filter == '.':
			return True
		else:
			return False
			
	def get_variant_type(self,caller,genotype):
		split_alt = self.alt.split(',')
		if genotype == '.' or genotype == './.' or genotype == '.|.':
			return 'uncalled_ambiguous'
		elif genotype == '0' or genotype == '0/0' or genotype == '0|0':
			return False
		else:
			alt = split_alt[int(genotype[-1:])-1]
			inequality_pattern = re.compile(r"^<\S+>$")
			if alt == '.':
				return False
			elif caller == 'PILON' and re.search(inequality_pattern,alt):
				return('structural')
			elif caller == 'GATK' and re.search(inequality_pattern,alt):
				return('inside_deletion')
			elif len(alt) == 1 and len(self.ref) == 1:
				return('SNP')
			elif len(alt) < len(self.ref):
				return('DELETION')
			elif len(alt) >= len(self.ref):
				return('INSERTION')
			else:
				return('unknown')
				
	def get_variant_length(self,genotype):
		if genotype == '.' or genotype == '0' or genotype == '0/0' or genotype == './.' or genotype == '0|0' or genotype == '.|.':
			return False
		else:
			ref_length = len(self.get_ref())
			genotype_end = re.search('(\d+)$',genotype)
			alt_length = len(self.get_alt(genotype_end.group(1)))
			if ref_length > 1 or alt_length > 1:
				return abs(alt_length - ref_length)
			else:
				return int(0)
			
# 	def get_genotype(self,index=0,min_gq=0,min_per_ad=float(0),min_tot_dp=0,return_flags=False):
# 		genotype = self.genotypes[index]
# 		parsed_genotype = genotype[0]   
# 		parsed_genotype_list = [parsed_genotype,0,0,0]
# 		### print(parsed_genotype) ####
# 		try:
# 			gq = self.get_GQ(parsed_genotype,index)
# 			if int(gq) < int(min_gq):
# 				parsed_genotype = '.'
# 				parsed_genotype_list[1] = 1
# 		except:
# 			pass
# 		try:
# 			percent_ad = self.get_percent_AD(index)
# 			if float(percent_ad) < min_per_ad:
# 				parsed_genotype = '.'
# 				parsed_genotype_list[2] = 1
# 		except:
# 			pass
# 		try:
# 			total_dp = self.get_total_DP(index)
# 			if int(total_dp) < int(min_tot_dp):
# 				parsed_genotype = '.'
# 				parsed_genotype_list[3] = 1
# 		except:
# 			pass
# 		if not return_flags:
# 			return parsed_genotype
# 		else:
# 			parsed_genotype_list[0] = parsed_genotype
# 			return parsed_genotype_list
			
	def get_genotype(self,index=0,min_gq=0,min_per_ad=float(0),min_tot_dp=0,het_binom_p=False,return_flags=False):  #### working on function to accomodate hets
		genotype = self.genotypes[index]
		parsed_genotype = genotype.split(':')[0]
		dip_flag = False
 		if len(parsed_genotype) == 3:
 			dip_flag = True
		parsed_genotype_list = [parsed_genotype,0,0,0,0]
		### print(parsed_genotype) ####
		try:
			gq = self.get_GQ(parsed_genotype,index)
			if int(gq) < int(min_gq):
				parsed_genotype = '.'
				if dip_flag == True:
					parsed_genotype = './.'
				parsed_genotype_list[1] = 1
		except:
			pass
		try:
			percent_ad = self.get_percent_AD(index)
			if float(percent_ad) < min_per_ad:
				parsed_genotype = '.'
				if dip_flag == True:
					parsed_genotype = './.'
				parsed_genotype_list[2] = 1
		except:
			pass
		try:
			total_dp = self.get_total_DP(index)
			if int(total_dp) < int(min_tot_dp):
				parsed_genotype = '.'
				if dip_flag == True:
					parsed_genotype = './.'
				parsed_genotype_list[3] = 1
		except:
			pass
		het_flag = self.is_het(index)
		if het_binom_p and het_flag:
			try:
				p = self.get_AD_binomial_p(index)
				if p < float(het_binom_p):
					parsed_genotype = './.'
					parsed_genotype_list[4] = 1
			except:
				pass
		if not return_flags:
			return parsed_genotype
		else:
			parsed_genotype_list[0] = parsed_genotype
			return parsed_genotype_list
	
	def is_het(self,index=0): #### currently works only on biallelic sites
		het = False
		genotype = self.genotypes[index]
		parsed_genotype = genotype.split(':')[0]
		gt_list = list(parsed_genotype)
		if '1' in gt_list and '0' in gt_list:
			het = True
		return het
		
	def get_GQ(self,parsed_genotype,index=0):
		gq = 'Undefined'
		fields = self.genotypes[index]
		gt_fields = fields.split(':')
		try:
			gq_index = self.get_GQ_index(parsed_genotype)
			gq = gt_fields[gq_index]
		except:
			pass
		return gq
		
	def get_percent_AD(self,index=0): #### currently works only on biallelic sites
		percent_AD = 'Undefined'
		fields = self.genotypes[index]
		gt_fields = fields.split(':')
		try:
			ad_index = self.get_AD_index()
			ad = gt_fields[ad_index]
			split_ad = ad.split(',')
			if split_ad[0] != '0':
				percent_AD = float(split_ad[1])/(float(split_ad[0])+float(split_ad[1]))
			elif split_ad[1] != '0':
				percent_AD = float(1)
		except:
			pass
		return percent_AD
		
	def get_AD_binomial_p(self,index=0): #### currently works only on biallelic sites
		pvalue = False
		fields = self.genotypes[index]
		gt_fields = fields.split(':')
		try:
			ad_index = self.get_AD_index()
			ad = gt_fields[ad_index]
			split_ad = ad.split(',')
			if split_ad[1] != '0':
				pvalue = binom_test(int(split_ad[1]), (int(split_ad[0]) + int(split_ad[1])))
			elif split_ad[0] != '0':
				pvalue = binom_test(int(split_ad[0]), (int(split_ad[1]) + int(split_ad[0])))
		except:
			pass
		return pvalue		
		
	def get_total_DP(self,index=0):
		total_DP = 'Undefined'
		fields = self.genotypes[index]
		## print(fields) ##
		gt_fields = fields.split(':')
		## print(gt_fields) ##
		try:
			dp_index = self.get_DP_index()
			total_DP = gt_fields[dp_index]
		except:
			pass
		return total_DP
		
	def get_GQ_index(self,parsed_genotype):
		gq_index = 'Undefined'
		fields = self.format
		format_fields = fields.split(':')
		### print(format_fields) ####
		### print(parsed_genotype) ####
		if parsed_genotype == '0' or parsed_genotype == '0/0' or parsed_genotype == '0|0':
			### print("try") ####
			try:
				gq_index = format_fields.index('RGQ')
			except:
				try:
					gq_index = format_fields.index('GQ')
				except:
					pass
		else:
			try:
				gq_index = format_fields.index('GQ')
			except:
				pass
		### print(gq_index) ###
		return gq_index
		
	def get_AD_index(self):
		ad_index = 'Undefined'
		fields = self.format
		format_fields = fields.split(':')
		try:
			ad_index = format_fields.index('AD')
		except:
			pass
		return ad_index
		
	def get_DP_index(self):
		ad_index = 'Undefined'
		fields = self.format
		format_fields = fields.split(':')
		try:
			ad_index = format_fields.index('DP')
		except:
			pass
		return ad_index	
	
	def get_chrom(self):
		return self.chrom

	def get_pos(self):
		return self.pos
		
	def get_id(self):
		return self.id

	def get_ref(self):
		return self.ref

	def get_alt(self,genotype):
		if genotype == '0' or genotype == '.' or genotype == '0/0' or genotype == './.' or genotype == '0|0' or genotype == '.|.':
			return False
		else:
			split_genotype = re.split(r"[/|]", genotype)
			if len(set(split_genotype)) > 1:
				return 'N'
			else:
				split_alt = self.alt.split(',')
				genotype_end = re.search('(\d+)$',genotype)
				index = int(genotype_end.group(1))-1
				return split_alt[index]		
			
	def get_alt_field(self):
		return self.alt
			
	def get_snpeff_annot(self, alt):
		fields = self.info.split(';')
		for field in fields:
			if re.match('ANN',field):
				ann_fields = field.split(',')
				## print ann_fields ##
				for ann_field in ann_fields:
					m = re.match('(ANN=)*(\w+)',ann_field)
					if m.group(2):
						if m.group(2) == alt:
							return ann_field.replace('ANN=','')				
		return False
		
	def get_snpeff_effect(self, snpeff_annot):
		try:
			annot_sections = snpeff_annot.split('|')
			return annot_sections[1]
		except:
			return False
		
	def get_snpeff_impact(self, snpeff_annot):
		try:
			annot_sections = snpeff_annot.split('|')
			return annot_sections[2]	
		except:
			return False
			
	def get_snpeff_feature(self, snpeff_annot):
		try:
			annot_sections = snpeff_annot.split('|')
			return annot_sections[3]	
		except:
			return False

	def get_qual(self):
		return self.qual

	def get_filter(self):
		return self.filter

	def get_info(self):
		return self.info
		
	def get_AF(self):
		AF = False
		fields = self.info.split(';')
		for field in fields:
			m = re.search('AF=([\d\.]+)',field)
			try:
				AF = float(m.group(1))
			except:
				pass
		return AF
		
	def get_QP(self):
		QP = False
		fields = self.info.split(';')
		for field in fields:
			m = re.search('QP=(\d+),(\d+),(\d+),(\d+)',field)
			try:
				QP = m.group(1,2,3,4)
				QP = [int(q) for q in QP]
			except:
				pass
		return QP
				
	def get_MAF_from_QP(self):
		MAF = False
		QP = self.get_QP()
		try:
			MAF = (100 - max(QP)) / 100
		except:
			pass
		return MAF

	def get_format(self):
		return self.format

	def get_genotypes_fields(self):
		return self.genotypes
		
	def get_genotypes_field(self,sample_index):
		return self.genotypes[sample_index]

	def get_vcf_annot(self):
		return self.vcf_annot
		
	def is_singleton(self):
		nonzeros = 0
		last_nonzero_index = False
		for i in range(0,len(self.genotypes)):
			current_genotype = self.get_genotype(index=i)
			if re.search(r"[1-9]",current_genotype):
				nonzeros += 1
				last_nonzero_index = i
		if nonzeros == 1:
			return str(last_nonzero_index)
		else:
			return False
			
	def is_biallelic(self):
		split_alt = self.alt.split(',')
		if len(split_alt) > 1:
			return False
		else:
			return True
			
	def count_ambig_genotypes(self):
		ambig = 0
		for current_genotype in self.genotypes:
			if current_genotype.split(':')[0] == '.' or current_genotype.split(':')[0] == './.':
				ambig += 1
		return ambig
		
	def get_genotype_profile(self):
		profile = list()
		for current_genotype in self.genotypes:
			profile.append(current_genotype.split(':')[0])
		return profile
		
class VcfHeader:
	"""
	
	"""
	def __init__(self, vcf_file):
		self.samples = []
		self.caller = False
		self.snpeff = False
		self.sample_columns = {}
		self.contigs = []

		comment_pattern = re.compile(r"^#")
		
		with open(vcf_file, 'r') as file:
			for full_line in file:
				line = full_line.rstrip()
				if (re.search(comment_pattern, line)):
					if (re.match('#CHROM', line)):
						fields = line.split('\t')
						for i in range(9,len(fields)):
							self.samples.append(fields[i])
							self.sample_columns[fields[i]] = i
					elif (re.match('##PILON', line)):
						self.caller = 'PILON'
					elif (re.match('##GATK', line)):
						self.caller = 'GATK'
					elif (re.match('##SnpEff', line)):
						self.snpeff = True
					elif (re.match('##contig',line)):
						m = re.search('##contig=<ID=([^,]+),',line)
						self.contigs.append(m.group(1))
				else:
					break
					
		if self.samples == ['SAMPLE']:
			self.samples = [vcf_file]
			self.sample_columns[vcf_file] = 9
			sys.stderr.write("No sample name in " + vcf_file + ", using file name.\n")
					
	def get_samples(self):
		return self.samples
		
	def get_caller(self):
		return self.caller
		
	def get_sample_index(self,sample):
		index = self.sample_columns[sample]-9
		return index
		
	def get_snpeff_status(self):
		return self.snpeff
		
	def get_contigs(self):
		return self.contigs
	


		
