from bs4 import BeautifulSoup
import json
import re
import os
from collections import defaultdict
from Bio import SeqIO, AlignIO
from Bio.Align.Applications import MafftCommandline
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

class GRec:
	""" data structures and methods to deal with individual records parsed from GISAID """
	def __init__(self,record_id,html):
		self.record_id = record_id
		self.html = html
		self.data = self.parse_html()

	def write_html_to_outfile(self):
		""" write the record to the given outfile """
		outfile = self.record_id + ".html"
		with open(outfile,'w') as o:
			o.write(self.html)

	def _parse_record_html_tr_tag_to_key_value_pair(self,tr_tag):
		""" take a tr tag from the record and return a key value pair if found otherwise None """
		tds = [x.text for x in list(tr_tag.children) if x != '\n']
		if (len(tds) != 2) or (':' not in tds[0]):
			return None
		tds[0] = tds[0].rstrip(':')
		return tds

	def parse_html(self):
		""" parse the html from a GISAID record into a data dict """
		soup = BeautifulSoup(self.html, 'html.parser')
		trs = soup.find_all('tr')
		data = {}
		address_count = 1
		for tr in trs:
			result = self._parse_record_html_tr_tag_to_key_value_pair(tr)
			if result is not None:
				## handle addresses differently, as there are multiple address fields
				if result[0] == 'Address':
					result[0] += str(address_count)
					address_count += 1	
				data[result[0]] = result[1]
			## add an empty "Comment" field if it wasn't found
			if 'Comment' not in data:
				data['Comment'] = ""
		fasta = soup.find_all('pre')[0].text.replace('\n','')
		fasta = re.sub(r'\s+$',r'',fasta)
		fasta = re.sub(r'^\s+',r'',fasta)
		fasta = re.sub(r'^(>.+EPI.+\d+)',r'\1\n',fasta)
		data['Fasta'] = fasta	
		return data

	def write_data_to_json_outfile(self):
		""" writes this object's data dict to a json outfile """
		outfile = self.record_id + '.json'
		with open(outfile,'w') as o:
			json.dump(self.data,o)

	def get_SeqRecord(self):
		""" take the fasta from the json and make a SeqIO object """
		fasta_header = self.data['Fasta'].split('\n')[0].lstrip(">")
		fasta_seq = self.data['Fasta'].split('\n')[1]
		epi_match = re.search(r'(.+)\|(EPI_\w+_\d+)',fasta_header)
		id = None
		if epi_match is not None:
			desc = epi_match.group(1)
			id = epi_match.group(2)
		return SeqRecord(Seq(fasta_seq),id=id, description=desc)

class GAlign:
	""" methods and data structures to perform alignments between a GRec fasta and the covid-19 
	genome sequence, and to produce vcf files from the alignment """

	def __init__(self,record_SeqRecord,reference_fasta_filepath="/home/aawan/SHARED/COVID-19/REF/NC_045512.fa"):
		self.record = record_SeqRecord
		self.reference = SeqIO.read(reference_fasta_filepath, "fasta") 
		self.reference_fasta_filepath = reference_fasta_filepath

	def _make_alignment_input_file(self):
		""" make an input file for the sequence alignment method based on the given record BioSeq object and the reference record """
		output_filepath = f'reference_plus_{self.record.id}.fa'
		SeqIO.write([self.reference, self.record], output_filepath, "fasta")
		return output_filepath

	def _align_fasta_to_ref(self,input_filepath):
		""" perform an alignment of the record fasta to the reference genome, write the alignment output, return the alignment object """
		mafft_cline = MafftCommandline(input=input_filepath)
		mafft_cline.clustalout = True
		stdout, stderr = mafft_cline()
		output_filepath = self.record.id + ".aln"
		with open(output_filepath,'w') as o:
			o.write(stdout)
		return output_filepath

	def process_record(self):
		""" wrapper method to call the above methods to align the fasta from the given record and to make a vcf file from the alignment """
		self._make_alignment_input_file()
		fasta_filepath = self._make_alignment_input_file()
		aln_filepath = self._align_fasta_to_ref(fasta_filepath)
		os.remove(fasta_filepath)
		a2v = AlnToVCF(self.reference_fasta_filepath,aln_filepath)
		vcf_filepath = self.record.id + ".vcf"
		a2v.write_vcf_from_mutations(vcf_filepath)


class GVCF:
	""" methods and data structures to represent and interact with VCF files made from GISAID records """

	def __init__(self,header,GVCFLine_list):
		self.header = header
		self.fields = ['CHROM','ID','POS','REF','ALT']
		self.GVCFLine_list = GVCFLine_list

	@classmethod
	def from_vcf_file(cls,vcf_filepath):
		""" create a GVCF object from a vcf file """
		with open(vcf_filepath) as v:
			header = v.readline().rstrip('\n')
			v.readline() # ignore the fields
			GVCFLine_list = []
			for line in v:
				GVCFLine_list.append(GVCFLine.from_line(line))
		return cls(header,GVCFLine_list)

	def write_vcf_file(self,outpath):
		""" write this object to a vcf file """
		with open(outpath,'w') as o:
			o.write(self.header + '\n')
			o.write('\t'.join(self.fields) + '\n')
			for gvcf_line in self.GVCFLine_list:
				o.write(gvcf_line.get_string_representation() + '\n')

	def get_line_list(self):
		""" return the lines of this object as a list of GVCF_Line objects """
		return self.GVCFLine_list


class GVCFLine:
	""" class to represent a single line (all properties should be strings) from a VCF file made from a GISAID record """
	
	def	__init__(self,chrom,id,pos,ref_seq,alt_seq):
		self.chrom = chrom
		self.id = id
		self.pos = pos
		self.ref_seq = ref_seq
		self.alt_seq = alt_seq

	@classmethod
	def from_line(cls,line):
		""" create a GVCFLine object from a line of a vcf file """
		tabs = line.rstrip('\n').split('\t')
		return cls(tabs[0],tabs[1],tabs[2],tabs[3],tabs[4])

	def get_string_representation(self):
		""" return a string representation of this object, suitable for printing or writing to a VCF file """
		return '\t'.join([self.chrom,self.id,str(self.pos),self.ref_seq,self.alt_seq])

	def get_int_position(self):
		""" return the integer representation of the position """
		return int(self.pos)

class AlnToVCF:
	""" methods and data structures to take an alignment file and a reference fasta file and produce a vcf file """
	
	def __init__(self,reference_fasta_filepath,alignment_filepath):
		alignment = AlignIO.read(open(alignment_filepath), "clustal")
		self.ref_seq_for_vcf = SeqIO.read(reference_fasta_filepath, "fasta").seq.upper()
		self.reference = alignment[0].seq.upper()
		self.record = alignment[1].seq.upper()
		self.mutations = self._get_mutations()
		self.mutations = self._consolidate_mutations(self.mutations)

	def write_vcf_from_mutations(self,vcf_filepath):
		""" write a vcf file using this object's mutations """
		header = "##fileformat=VCFv4.2"
		GVCFLine_list = []
		CHROM = "NC_045512"
		ID = "."
		for pos,mut in self.mutations.items(): # keep in mind that the self.mutations objects are zero-indexed, whereas VCF is 1-indexed
			if mut.type == 'SUB':
				POS = mut.ref_start + 1
				REF = mut.ref_seq
				ALT = mut.alt_seq		
			elif mut.type == 'UNK':
				POS = mut.ref_start + 1
				REF = mut.ref_seq
				ALT = mut.alt_seq			
			elif mut.type == 'DEL':
				POS = mut.ref_start
				ALT = self.ref_seq_for_vcf[POS:POS+1]
				REF = ALT + mut.ref_seq
			elif mut.type == 'INS':
				POS = mut.ref_start
				REF = self.ref_seq_for_vcf[POS:POS+1]
				ALT = REF + mut.alt_seq
			GVCFLine_list.append( GVCFLine(CHROM,ID,str(POS),str(REF),str(ALT)) )
		gvcf = GVCF(header,GVCFLine_list)
		gvcf.write_vcf_file(vcf_filepath)
			
	def _get_mutations(self):
		""" iterate through the reference and record sequences and return a list of AlnMut objects """
		muts = defaultdict(list)
		ref_pos = 0
		for i in range(len(self.reference) - 1):
			ref_seq = self.reference[i:i+1]
			rec_seq = self.record[i:i+1]
			if ref_seq == rec_seq:
				ref_pos += 1
				continue
			type = None
			## determine the mutation type, add an AlnMut object to the muts list, and update positions accordingly
			current_mut = None
			if (ref_seq != '-') and (rec_seq != '-'): # SUB or UNK
				if rec_seq == 'N':
					type = 'UNK'
				else:
					type = 'SUB'
				current_mut = AlnMut(type,ref_seq,rec_seq,ref_pos) 
			elif rec_seq == '-': # DEL
				type = 'DEL'
				current_mut = AlnMut(type,ref_seq,rec_seq,ref_pos)
			elif ref_seq == '-': # INS
				type = 'INS'
				current_mut = AlnMut(type,ref_seq,rec_seq,ref_pos)
			#print(ref_pos,current_mut)
			muts[ref_pos].append(current_mut)
			if type != 'INS':
				ref_pos += 1
		return muts

	def _consolidate_mutations(self,mutations):
		""" given a list of AlnMut object and return a list of consolidated AlnMut objects """
		consolidated = defaultdict(AlnMut)
		## first consolidate the INS
		for pos in sorted(mutations):
			if len(mutations[pos]) > 1:
				ins_ref_start = mutations[pos][0].ref_start
				ins_alt_seq = ''.join([str(x.alt_seq) for x in mutations[pos]])
				ins_ref_seq = '-'
				consolidated[pos] = AlnMut('INS',ins_ref_seq,ins_alt_seq,ins_ref_start)
			else:
				consolidated[pos] = mutations[pos][0]
		## now consolidate the DEL and UNK
		consolidation = False	
		first_time = True
		while consolidation or first_time:
			first_time = False
			consolidation = False
			sorted_pos = list(sorted(consolidated))
			for i in range(len(sorted_pos) - 1):
				cons_mut = AlnMut.consolidate(consolidated[sorted_pos[i]], consolidated[sorted_pos[i+1]])
				#print(len(sorted_pos),i,sorted_pos[i],cons_mut)
				if cons_mut is None:
					continue
				#print(cons_mut)
				consolidated[sorted_pos[i]] = cons_mut
				del(consolidated[sorted_pos[i+1]])
				consolidation = True
				break
		return consolidated				 
				

class AlnMut:
	""" data structure to represent a mutation of an alignment of a record relative to a reference sequence """
	
	def __init__(self,mut_type,ref_seq,alt_seq,ref_start):
		self.type = mut_type # one of SUB,INS,DEL,UNK describing the relationship of the record seq to the reference seq for this mutation
		self.ref_seq = ref_seq # the sequence of the reference between ref_start and ref_end positions
		self.alt_seq = alt_seq # the sequence of the aligned record between ref_start and ref_end positions
		self.ref_start = ref_start # the starting position on the reference of this mutation
		## calculate the position of the end of the mutation with regard to the reference
		self.ref_end = self.ref_start + len(self.ref_seq)
		if self.type == 'INS':
			self.ref_end = self.ref_start + 1

	def __str__(self):
		""" custom print function """
		return ",".join([self.type,str(self.ref_seq),str(self.alt_seq),str(self.ref_start),str(self.ref_end)])

	@classmethod
	def consolidate(cls,mut1,mut2):
		""" constructor to produce a consolidated mutation from two compatible,contiguous mutations """
		if (mut1.type != mut2.type) or (mut1.type == 'SUB'):
			return None
		if mut1.ref_end != mut2.ref_start:
			return None
		ret_type = mut1.type
		ret_ref_start = mut1.ref_start
		ret_ref_seq = mut1.ref_seq
		ret_alt_seq = mut1.alt_seq
		if (ret_type == 'DEL') or (ret_type == 'UNK'):
			ret_ref_seq += mut2.ref_seq
		if ret_type == 'UNK':
			ret_alt_seq += mut2.alt_seq			
		return cls(ret_type, ret_ref_seq, ret_alt_seq, ret_ref_start)

if __name__ == "__main__":
	pass
