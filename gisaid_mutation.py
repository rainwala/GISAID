from gisaid_record import GVCFLine
from Bio import SeqIO

class GisVar:
	""" methods and data structures to process and represent variants from vcf files made from alignments of records to the reference 
	with nomenclature taken from the Human Genome Variation Society (HGVS) standards 
	(https://onlinelibrary.wiley.com/doi/pdf/10.1002/humu.22981)"""
	
	def __init__(self,GVCFLine,record_id_list):
		self.HGVS_types = {
										'substitution':'>', 
										'insertion':'ins',
										'deletion':'del',
		}
		self.sc2p = SARSCov2Prot()
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

	def _construct_genomic_variant_name(self):
		""" make the name of this variant with respect to the reference genome """
		name = None
		if self.type == 'substitution':
			name = f'g.{self.genome_start}{self.GVCFLine.ref_seq}{self.HGVS_types[self.type]}{self.GVCFLine.alt_seq}'	
		elif  self.type == 'deletion':
			name = f'g.{self.genome_start}_{self.genome_end}{self.HGVS_types[self.type]}'
		elif self.type == 'insertion':
			## the '[1:]' index on alt_seq is because we need to trim the first base from the vcf alt seq
			name = f'g.{self.genome_start}_{self.genome_end}{self.HGVS_types[self.type]}{self.GVCFLine.alt_seq[1:]}'
		return name		

	def _construct_protein_variant_name(self):
		""" make the name of this variant with respect to a viral protein, if the variant is in a protein """
		return self.sc2p.get_var_protein_bounds_from_genomics_bounds(self.genome_start,self.genome_end)

	def _determine_synonymity(self):
		""" based on the protein name, determine whether this mutation is synonymous or not """
		pass


class SARSCov2Prot:
	""" methods and data strcutures to represent and interact with proteins from the SARS-Cov2 genome """

	def __init__(self,reference_fasta_filepath="/home/aawan/SHARED/COVID-19/REF/NC_045512.fa"):
		self.reference = SeqIO.read(reference_fasta_filepath, "fasta")
		self.reference_protein_locations = {
			'nsp1':[266,805], 'nsp2':[806,2719], 'nsp3':[2720,8554], 'nsp4':[8555,10054], '3C-like proteinase':[10055,10972],
			'nsp6':[10973,11842], 'nsp7':[11843,12091], 'nsp8':[12092,12685], 'nsp9':[12686,13024], 'nsp10':[13025,13441],
			'nsp11':[13442,13480], 'nsp12':[13442,16236], 'nsp13':[16237,18039], 'nsp14': [18040,19620], 'nsp15': [19621,20658],
			'nsp16':[20659,21552], 'Spike protein':[21563,25384], 'ORF3a':[25393,26220], 'E':[26245,26472], 'M':[26523,27191], 'ORF6':[27202,27387],
			'ORF7a':[27394,27759], 'ORF7b':[27756,27887], 'ORF8':[27894,28259], 'N':[28274,29533], 'ORF10':[29558,29674],
		}

	def get_var_protein_bounds_from_genomics_bounds(self,genome_start,genome_end):
		""" given the starting and ending genomic position of a variant, 
		return the protein (if any) that span occurs in, and the  protein co-ords """
		start_prot = None
		start_prot_start = None
		end_prot = None
		end_prot_end = None
		for prot,coords in self.reference_protein_locations.items():
			if (genome_start >= coords[0]) and (genome_start <= coords[1]):
				start_prot = prot
				start_prot_start = (genome_start - coords[0]) + 1
			if (genome_end >= coords[0]) and (genome_end <= coords[1]):
				end_prot = prot
				end_prot_end = (genome_end - coords[0]) + 1
		return [start_prot,start_prot_start,end_prot,end_prot_end]
