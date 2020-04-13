from gisaid_record import GVCFLine
from Bio import SeqIO
from Bio.SeqUtils import seq3
import math

class GisVar:
	""" methods and data structures to process and represent variants from vcf files made from alignments of records to the reference 
	with nomenclature taken from the Human Genome Variation Society (HGVS) standards 
	(https://onlinelibrary.wiley.com/doi/pdf/10.1002/humu.22981)"""
	
	def __init__(self,GVCFLine,record_id_list):
		self.reference_genome_name = 'NC_045512.2'
		self.reference_fasta_filepath = "/home/aawan/SHARED/COVID-19/REF/NC_045512.fa"
		self.reference = SeqIO.read(self.reference_fasta_filepath, "fasta")
		self.reference_protein_locations = {
      'nsp1':[266,805], 'nsp2':[806,2719], 'nsp3':[2720,8554], 'nsp4':[8555,10054], '3C-like proteinase':[10055,10972],
      'nsp6':[10973,11842], 'nsp7':[11843,12091], 'nsp8':[12092,12685], 'nsp9':[12686,13024], 'nsp10':[13025,13441],
      'nsp11':[13442,13480], 'nsp12':[13442,16236], 'nsp13':[16237,18039], 'nsp14': [18040,19620], 'nsp15': [19621,20658],
      'nsp16':[20659,21552], 'Spike protein':[21563,25384], 'ORF3a':[25393,26220], 'E':[26245,26472], 'M':[26523,27191], 'ORF6':[27202,27387],
      'ORF7a':[27394,27759], 'ORF7b':[27756,27887], 'ORF8':[27894,28259], 'N':[28274,29533], 'ORF10':[29558,29674],
		}
		self.HGVS_types = {
										'substitution':'>', 
										'insertion':'ins',
										'deletion':'del',
		}
		self.GVCFLine = GVCFLine
		self.type = self._get_type()
		self.genome_start, self.genome_end = self._get_genomic_bounds()
		self.record_ids = record_id_list

	def _get_type(self):
		""" return the type of variant """
		if len(self.GVCFLine.ref_seq) == len(self.GVCFLine.alt_seq):
			return 'substitution'
		if (len(self.GVCFLine.alt_seq) == 1) and (len(self.GVCFLine.ref_seq) > len(self.GVCFLine.alt_seq)):
			return 'deletion'
		if (len(self.GVCFLine.ref_seq) == 1) and (len(self.GVCFLine.alt_seq) > len(self.GVCFLine.ref_seq)):
			return 'insertion'
	
	def _get_genomic_bounds(self):
		""" return the genome start and end position """
		return self.GVCFLine.get_int_position(), self.GVCFLine.get_int_position() + len(self.GVCFLine.ref_seq)

	def make_genomic_variant_name(self):
		""" return an hgvs object representing this variant with respect to the reference genome """
		name = None
		if self.type == 'substitution':
			name = f'{self.reference_genome_name}:g.{self.genome_start}{self.GVCFLine.ref_seq}{self.HGVS_types[self.type]}{self.GVCFLine.alt_seq}'	
		elif  self.type == 'deletion':
			name = f'{self.reference_genome_name}:g.{self.genome_start}_{self.genome_end}{self.HGVS_types[self.type]}'
		elif self.type == 'insertion':
			## the '[1:]' index on alt_seq is because we need to trim the first base from the vcf alt seq
			name = f'{self.reference_genome_name}:g.{self.genome_start}_{self.genome_end}{self.HGVS_types[self.type]}{self.GVCFLine.alt_seq[1:]}'
		return name

	def make_protein_variant_name(self):
		""" make the name of this variant with respect to a viral protein, if the variant is in a protein """
		overlaps = self._get_proteins_overlapping_variant()
		if len(overlaps) != 1:
			return None
		prot_name = next(iter(overlaps))
		ref_AA_seq = self._get_protein_affected_reference_sequence(prot_name)
		## handle the easy case where this is a frameshift
		ref_first_AA_3letter = seq3(ref_AA_seq[0])
		var_start_codon, var_end_codon = self._get_AA_coords_within_protein(prot_name)
		if (self.type != 'substitution'):
			indel_length = abs(len(self.GVCFLine.alt_seq) - len(self.GVCFLine.ref_seq))
			if indel_length % 3 != 0:
				return 'p.(' + ref_first_AA_3letter + str(var_start_codon) + 'fs)'
		ref_last_AA_3letter = seq3(ref_AA_seq[-1])	
		name = None
		## handle substitutions
		if self.type == 'substitution':
			var_AA_seq = self._get_protein_variant_sequence(prot_name)
			name =  'p.' + ref_first_AA_3letter + str(var_start_codon) + seq3(var_AA_seq[0])
		## handle deletions
		elif self.type == 'deletion':
			name = 'p.' + ref_first_AA_3letter + str(var_start_codon) + "_" + ref_last_AA_3letter + str(var_end_codon) + 'del'
		## handle insertions
		elif self.type == 'insertion':
			var_AA_seq = self._get_protein_variant_sequence(prot_name)
			name = 'p.' + ref_first_AA_3letter + str(var_start_codon) + "_" + ref_last_AA_3letter + str(var_end_codon) + 'ins' + seq3(var_AA_seq)
		return name

	def is_synonymous(self):
		""" based on the protein name, return whether this mutation is synonymous or not """
		pass

	def _get_proteins_overlapping_variant(self):
		""" return the protein(s) (if any) that the variant occurs in """
		overlaps = set()
		for prot,coords in self.reference_protein_locations.items():
			if (self.genome_start >= coords[0]) and (self.genome_start <= coords[1]):
				overlaps.add(prot)
			if (self.genome_end >= coords[0]) and (self.genome_end <= coords[1]):
				overlaps.add(prot)
		return overlaps

	def _get_DNA_coords_within_protein(self,protein_name):
		""" return the start and end DNA position of the variant with respect to the protein start position """
		prot_start = self.reference_protein_locations[protein_name][0]
		variant_DNA_start_in_protein = (self.genome_start - prot_start) + 1
		variant_DNA_end_in_protein = (self.genome_end - prot_start) + 1
		return variant_DNA_start_in_protein,variant_DNA_end_in_protein

	def _get_AA_coords_within_protein(self,protein_name):
		""" return the start and end AA position of the variant with respect to the protein start position """
		variant_DNA_start_in_protein,variant_DNA_end_in_protein = self._get_DNA_coords_within_protein(protein_name)	
		return math.ceil(variant_DNA_start_in_protein/3), math.ceil(variant_DNA_end_in_protein/3)

	def _convert_protein_DNA_coords_to_genomic(self,protein_name,prot_DNA_start,prot_DNA_end):
		""" return the genomic coords given DNA coords with respect to the protein start position"""
		genome_start_pos = self.reference_protein_locations[protein_name][0] + prot_DNA_start
		genome_end_pos = self.reference_protein_locations[protein_name][0] + prot_DNA_end
		return genome_start_pos,genome_end_pos

	def _get_protein_affected_reference_sequence(self,protein_name,verbose=False):
		""" given the protein name, return the translated reference sequence, 
		starting from the entire codon containing 'self.genome_start' and ending with 
		the entire codon containing 'self.genome_end' """
		prot_DNA_start, prot_DNA_end = self._get_DNA_coords_within_protein(protein_name)
		## try to extend the start and end co-ords, depending on the frame
		start_adjustment = prot_DNA_start % 3
		end_adjustment = prot_DNA_end % 3
		ref_prot_DNA_start = prot_DNA_start - start_adjustment
		ref_prot_DNA_end = prot_DNA_end
		if end_adjustment > 0:
			ref_prot_DNA_end += (3 - end_adjustment)
		## convert the co-ords back to genomic ones
		adjusted_genome_start,adjusted_genome_end = self._convert_protein_DNA_coords_to_genomic(protein_name,ref_prot_DNA_start,ref_prot_DNA_end)
		## return the translated sequence, making sure to zero-index the sequence that is retreived
		## we need to add the following to allow the variant name in hgvs format to include the next amino acid in the reference
		if self.type == 'insertion':
			adjusted_genome_end += 3
		ref_DNA_seq = self.reference.seq[adjusted_genome_start - 1 : adjusted_genome_end - 1]
		if not verbose:
			return ref_DNA_seq.translate()
		return ref_DNA_seq, ref_prot_DNA_start, ref_prot_DNA_end, prot_DNA_start, prot_DNA_end

	def _get_protein_variant_sequence(self,protein_name):
		""" given the protein name, return the translated variant sequence, 
    starting from the entire codon containing 'self.genome_start' and ending with 
    the entire codon containing 'self.genome_end' """
		## get the reference DNA sequence and start and end DNA positions of the affected reference codons
		ref_DNA_seq, ref_prot_DNA_start, ref_prot_DNA_end, prot_DNA_start, prot_DNA_end = self._get_protein_affected_reference_sequence(protein_name,verbose=True)
		print('ref',ref_DNA_seq)
		if self.type == 'insertion':
			ins_pos = prot_DNA_start - ref_prot_DNA_start
			var_DNA_seq = ref_DNA_seq[0:ins_pos] + self.GVCFLine.alt_seq + ref_DNA_seq[ins_pos + 1:]
		elif self.type == 'substitution':
			sub_pos = prot_DNA_start - ref_prot_DNA_start # sub_pos should only ever be 0,1, or 2
			print('sub',sub_pos,ref_DNA_seq[0:sub_pos - 1],self.GVCFLine.alt_seq,ref_DNA_seq[sub_pos + 1:])
			var_DNA_seq = ref_DNA_seq[0:sub_pos] + self.GVCFLine.alt_seq
			if sub_pos < 2:
				var_DNA_seq += ref_DNA_seq[sub_pos + 1:]
		print(self.type,var_DNA_seq)
		return var_DNA_seq.translate()
			
