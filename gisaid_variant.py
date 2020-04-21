from gisaid_record import GVCFLine
from Bio import Seq,SeqIO
from Bio.SeqUtils import seq3
import math
import difflib
import re

class GisVar:
	""" parent class for the GisSub, GisDel, and GisIns classes below: 
	methods and data structures to process and represent variants from vcf files made from alignments of records to the reference 
  with nomenclature taken from the Human Genome Variation Society (HGVS) standards 
  (https://onlinelibrary.wiley.com/doi/pdf/10.1002/humu.22981) """
	
	def __init__(self,GVCFLine,record_id_list):
		self.reference_genome_name = 'NC_045512.2'
		self.reference_fasta_filepath = "/home/aawan/SHARED/COVID-19/REF/NC_045512.fa"
		self.reference = SeqIO.read(self.reference_fasta_filepath, "fasta")
		self.reference_protein_locations = {
      'nsp1':[266,805], 'nsp2':[806,2719], 'nsp3':[2720,8554], 'nsp4':[8555,10054], '3C-like proteinase':[10055,10972],
      'nsp6':[10973,11842], 'nsp7':[11843,12091], 'nsp8':[12092,12685], 'nsp9':[12686,13024], 'nsp10':[13025,13441],
      'nsp11':[13442,13480], 'nsp12':[13442,16236], 'nsp13':[16237,18039], 'nsp14': [18040,19620], 'nsp15': [19621,20658],
      'nsp16':[20659,21552], 'S':[21563,25384], 'ORF3a':[25393,26220], 'E':[26245,26472], 'M':[26523,27191], 'ORF6':[27202,27387],
      'ORF7a':[27394,27759], 'ORF7b':[27756,27887], 'ORF8':[27894,28259], 'N':[28274,29533], 'ORF10':[29558,29674],
		}		 
		self.GVCFLine = GVCFLine
		self.record_id_list = record_id_list
		self.type = self.get_type()
		self.type_object_dict = {
			'sub':GisSub, 'del':GisDel, 'ins':GisIns, 
		}

	def get_type(self):
		""" return the type of variant this is """
		if len(self.GVCFLine.ref_seq) == len(self.GVCFLine.alt_seq):
			return 'sub'
		elif (len(self.GVCFLine.alt_seq) == 1) and (len(self.GVCFLine.ref_seq) > len(self.GVCFLine.alt_seq)):
			return 'del'
		elif (len(self.GVCFLine.ref_seq) == 1) and (len(self.GVCFLine.alt_seq) > len(self.GVCFLine.ref_seq)):
			return 'ins'
		
	def get_var_type_child(self):
		""" return the a child object of the type of variant this is """
		return self.type_object_dict[self.type](self.GVCFLine,self.record_id_list)

	def make_genomic_variant_name(self):
		""" make the name of this variant with respect to the reference genome """
		pass	

	def get_protein_name(self):
		""" return the name of the protein this variant is found in, or None """
		start_protein_name = self._get_protein_name_for_genomic_position(self.GVCFLine.pos)
		if start_protein_name is None:
			return None
		end_protein_name = self._get_protein_name_for_genomic_position(
			self.GVCFLine.pos + max(len(self.GVCFLine.ref_seq), len(self.GVCFLine.alt_seq)) 
		)
		if end_protein_name != start_protein_name:
			return None
		return start_protein_name

	def make_protein_variant_name(self):
		""" make the name of this variant with respect to a viral protein, if the variant is in a protein """
		pass
	
	def is_synonymous(self,protein_variant_name):
		""" based on the protein-level variant name, return whether this mutation is synonymous or not """
		if protein_variant_name is None:
			return True
		if ('fs' in protein_variant_name) or ('del' in protein_variant_name) or ('ins' in protein_variant_name):
			return False
		sub_match = re.search(r':p\.([A-Z][a-z]{2})(\d+)([A-Z][a-z]{2})',protein_variant_name)
		if sub_match is not None:
			if sub_match.group(1) == sub_match.group(3):
				return True
			else:
				return False
		return None

	def _get_protein_name_for_genomic_position(self,genome_position):
		""" return the name of the protein this genome position falls in, or None """
		for prot,coords in self.reference_protein_locations.items():
			if (genome_position >= coords[0]) and (genome_position <= coords[1]):
				return prot
		return None

	def _transform_genomic_position_to_protein(self,genome_position):
		""" get the 1-indexed within-protein position for the given 1-indexed genomic position """
		protein_name = self._get_protein_name_for_genomic_position(genome_position)
		if (protein_name is None) or (protein_name not in self.reference_protein_locations):
			return None
		return (genome_position - self.reference_protein_locations[protein_name][0]) + 1

	def _get_codon_number_for_genomic_position(self,genome_position):
		""" get the within protein codon number given a 1-indexed genomic position """
		protein_name = self._get_protein_name_for_genomic_position(genome_position)
		if (protein_name is None) or (protein_name not in self.reference_protein_locations):
			return None
		return math.ceil( self._transform_genomic_position_to_protein(genome_position) / 3 )

	def _get_frame_for_genomic_position(self,genome_position):
		""" get the within protein frame for a 1-indexed genomic position """
		protein_name = self._get_protein_name_for_genomic_position(genome_position)
		if (protein_name is None) or (protein_name not in self.reference_protein_locations):
			return None
		frame = self._transform_genomic_position_to_protein(genome_position) % 3
		if frame == 0:
			frame = 3
		return frame

	def _get_protein_codon_DNA_seq(self,protein_name,codon_number):
		""" given a protein name and a codon nummber, return the DNA sequence for that codon """
		if (protein_name is None) or (protein_name not in self.reference_protein_locations):
			return None
		prot_start_pos = self.reference_protein_locations[protein_name][0]
		codon_start_pos = prot_start_pos + 3*(codon_number - 1)
		codon_end_pos = codon_start_pos + 2
		return self.reference.seq[codon_start_pos - 1 : codon_end_pos].upper()
		
	def _get_ref_DNA_full_codon_seqs_between_two_genome_positions(self,genome_position1,genome_position2):
		""" get the reference AA sequence with which a comparison for naming will be made """
		codon_number1 = self._get_codon_number_for_genomic_position(genome_position1)
		codon_number2 = self._get_codon_number_for_genomic_position(genome_position2)
		## return None if either genomic position is outside a protein, or they are in different proteins
		if (codon_number1 is None) or (codon_number2 is None): 
			return None
		protein_name = self._get_protein_name_for_genomic_position(genome_position1)
		if protein_name != self._get_protein_name_for_genomic_position(genome_position2):
			return None
		ref_DNA_seq = ""
		for cn in range(codon_number1,codon_number2 + 1):
			ref_DNA_seq += self._get_protein_codon_DNA_seq(protein_name,cn)
		return ref_DNA_seq

	def _get_ref_and_var_AA_difference_for_indel_name(self,genomic_indel_type,start_codon_number,ref_AA_seq,var_AA_seq):
		""" compare the ref and var AA seqs to allow protein-level naming for indels """
		insertions = {}
		deletions = {}		
		## the offset is needed for cases where there is a protein delins resulting from a genomic deletion that deletes the first reference codon
		offset = 0
		if not var_AA_seq.startswith(ref_AA_seq[0]) and genomic_indel_type == 'del':
			offset = 1
		for i,s in enumerate(difflib.ndiff(ref_AA_seq, var_AA_seq)):
			if s[0]==' ': continue
			if s[0]=='-':
				deletions[(start_codon_number + i) - offset] = seq3(s[-1])
			elif s[0]=='+':
				insertions[(start_codon_number + i) - offset] = seq3(s[-1])
		## process the deletions
		name = ""
		if (len(deletions) > 0):
			if len(deletions) == 1:
				pos = sorted(deletions)[0]
				name = deletions[pos] + str(pos)
			else:
				min_pos = sorted(deletions)[0]
				max_pos = sorted(deletions)[-1]
				name = deletions[min_pos] + str(min_pos) + "_" + deletions[max_pos] + str(max_pos) 
			name += 'del'
		elif genomic_indel_type == 'ins':
			name = seq3(ref_AA_seq[0]) + str(start_codon_number) + '_' + seq3(ref_AA_seq[1]) + str(start_codon_number + 1)
		if len(insertions) > 0:
			name += 'ins'
			for pos in sorted(insertions):
				name += insertions[pos]
		return name


class GisSub(GisVar):
	""" methods and data structures to process and represent substitution variants """
	
	def __init__(self,GVCFLine,record_id_list):
		super().__init__(GVCFLine,record_id_list)
		## for a substitution, the start and end position are the same. These positions are 1-indexed
		self.genome_start = self.GVCFLine.get_int_position()
		self.genome_end = self.GVCFLine.get_int_position()

	def make_genomic_variant_name(self):
		""" make the name of this variant with respect to the reference genome """
		return f'{self.reference_genome_name}:g.{self.genome_start}{self.GVCFLine.ref_seq}>{self.GVCFLine.alt_seq}'

	def make_protein_variant_name(self):
		""" make the name of this variant with respect to a viral protein, if the variant is in a protein """
		protein_name = self._get_protein_name_for_genomic_position(self.genome_start)
		if protein_name is None:
			return None
		codon_number = self._get_codon_number_for_genomic_position(self.genome_start)
		ref_codon_DNA = self._get_protein_codon_DNA_seq(protein_name,codon_number)
		ref_codon_AA_3letter = seq3(ref_codon_DNA.translate())
		var_frame = self._get_frame_for_genomic_position(self.genome_start) 
		var_codon_AA_3letter = seq3( (ref_codon_DNA[:var_frame-1] + self.GVCFLine.alt_seq + ref_codon_DNA[var_frame:]).translate() )
		return f'{protein_name}:p.{ref_codon_AA_3letter}{codon_number}{var_codon_AA_3letter}'
		

class GisDel(GisVar):
	""" methods and data structures to process and represent deletion variants """
	
	def __init__(self,GVCFLine,record_id_list):
		super().__init__(GVCFLine,record_id_list)
		## for a deletion, the start postion is one after the position in the VCF file. These positions are 1-indexed
		self.genome_start = self.GVCFLine.get_int_position() + 1
		self.genome_end = self.genome_start  + (len(self.GVCFLine.ref_seq) - 2)

	def make_genomic_variant_name(self):
		""" make the name of this variant with respect to the reference genome """
		if self.genome_start == self.genome_end: # handle the case where the deletion is just 1 bp
			return f'{self.reference_genome_name}:g.{self.genome_start}del'
		return f'{self.reference_genome_name}:g.{self.genome_start}_{self.genome_end}del'
		
	def make_protein_variant_name(self):
		""" make the name of this variant with respect to a viral protein, if the variant is in a protein """
		protein_name = self._get_protein_name_for_genomic_position(self.genome_start)
		if protein_name is None:
			return None
		## if the variant spance 2 proteins, return None
		protein_name2 = self._get_protein_name_for_genomic_position(self.genome_end)
		if protein_name != protein_name2:
			return None
		start_codon_number = self._get_codon_number_for_genomic_position(self.genome_start)
		start_codon_DNA = self._get_protein_codon_DNA_seq(protein_name,start_codon_number)
		start_codon_AA_3letter = seq3(start_codon_DNA.translate())
		## if the variant is in a coding region and its length is not divisible by 3, return a frameshift name
		if (len(self.GVCFLine.ref_seq) - 1) % 3 != 0: # keep in mind the ref seq from a VCF always has one extra 3' base
			return f'{protein_name}:p.({start_codon_AA_3letter}{start_codon_number}fs)'
		ref_DNA_seq = self._get_ref_DNA_full_codon_seqs_between_two_genome_positions(self.genome_start,self.genome_end)
		ref_AA_seq = ref_DNA_seq.translate()
		start_frame = self._get_frame_for_genomic_position(self.genome_start)
		end_frame = self._get_frame_for_genomic_position(self.genome_end)
		del_start = start_frame-1
		num_codons_in_ref_DNA_Seq = (len(ref_DNA_seq))/3
		del_end = 3*int(num_codons_in_ref_DNA_Seq - 1) + end_frame
		var_DNA_seq = ref_DNA_seq[:del_start] + ref_DNA_seq[del_end:]
		var_AA_seq = var_DNA_seq.translate()
		indel_name = self._get_ref_and_var_AA_difference_for_indel_name('del',start_codon_number,ref_AA_seq,var_AA_seq)
		return f'{protein_name}:p.{indel_name}'


class GisIns(GisVar):
	""" methods and data structures to process and represent insertion variants """

	def __init__(self,GVCFLine,record_id_list):
		super().__init__(GVCFLine,record_id_list)
		## for an insertion, the start postion is one after the position in the VCF file. These positions are 1-indexed
		self.genome_start = self.GVCFLine.get_int_position()
		self.genome_end = self.genome_start + 1

	def make_genomic_variant_name(self):
		""" make the name of this variant with respect to the reference genome """
		return f'{self.reference_genome_name}:g.{self.genome_start}_{self.genome_end}ins{self.GVCFLine.alt_seq[1:]}'

	def make_protein_variant_name(self):
		""" make the name of this variant with respect to a viral protein, if the variant is in a protein """
		protein_name = self._get_protein_name_for_genomic_position(self.genome_start)
		if protein_name is None:
			return None
		start_codon_number = self._get_codon_number_for_genomic_position(self.genome_start)
		start_codon_DNA = self._get_protein_codon_DNA_seq(protein_name,start_codon_number)
		start_codon_AA_3letter = seq3(start_codon_DNA.translate())
		## if the variant is in a coding region and its length is not divisible by 3, return a frameshift name
		if (len(self.GVCFLine.alt_seq) - 1) % 3 != 0: # keep in mind the ref seq from a VCF always has one extra 3' base
			return f'{protein_name}:p.({start_codon_AA_3letter}{start_codon_number}fs)'
		ref_DNA_seq = self._get_ref_DNA_full_codon_seqs_between_two_genome_positions(self.genome_start,self.genome_end + 3)
		ref_AA_seq = ref_DNA_seq.translate()
		start_frame = self._get_frame_for_genomic_position(self.genome_start)
		ins_start = start_frame-1
		var_DNA_seq = ref_DNA_seq[:ins_start] + self.GVCFLine.alt_seq + ref_DNA_seq[ins_start+1:]
		var_AA_seq = var_DNA_seq.translate()
		indel_name = self._get_ref_and_var_AA_difference_for_indel_name('ins',start_codon_number,ref_AA_seq,var_AA_seq)
		return f'{protein_name}:p.{indel_name}'
