from gisaid_variant import GisVar,GisSub,GisDel,GisIns
from gisaid_record import GVCFLine
import pytest

class TestGisVar:

	@pytest.fixture
	def sub_gvcf(self):
		sub_gvcf_line = GVCFLine.from_line("NC_045512\t.\t21589\tC\tT\n")
		return GisVar(sub_gvcf_line,['EPI_ISL_421184'])

	@pytest.fixture
	def del_gvcf(self):
		del_gvcf_line = GVCFLine.from_line("NC_045512\t.\t21989\tGTTT\tG\n")
		return GisVar(del_gvcf_line,['EPI_ISL_421184'])

	@pytest.fixture
	def ins_gvcf(self):
		ins_gvcf_line = GVCFLine.from_line("NC_045512\t.\t3336\tA\tAATC\n")
		return GisVar(ins_gvcf_line,['EPI_ISL_410721'])

	def test_get_var_type_child(self,sub_gvcf,del_gvcf,ins_gvcf):
		assert GisSub == type(sub_gvcf.get_var_type_child())
		assert GisDel == type(del_gvcf.get_var_type_child())
		assert GisIns == type(ins_gvcf.get_var_type_child())
	
	@pytest.mark.parametrize("genome_position,expected",[(200,None),(266,'nsp1'),(29674,'ORF10')])
	def test_get_protein_name_for_genomic_position(self,sub_gvcf,genome_position,expected):
		assert expected == sub_gvcf._get_protein_name_for_genomic_position(genome_position)

	@pytest.mark.parametrize("genome_position,expected",[(200,None),(266,1),(29674,117)])
	def test_transform_genomic_position_to_protein(self,sub_gvcf,genome_position,expected):
		assert expected == sub_gvcf._transform_genomic_position_to_protein(genome_position)
	
	@pytest.mark.parametrize("genome_position,expected",[(200,None),(266,1),(16241,2),(29674,39)])
	def test_get_codon_number_for_genomic_position(self,sub_gvcf,genome_position,expected):
		assert expected == sub_gvcf._get_codon_number_for_genomic_position(genome_position)

	@pytest.mark.parametrize("genome_position,expected",[(200,None),(266,1),(16241,2),(29674,3)])
	def test_get_frame_for_genomic_position(self,sub_gvcf,genome_position,expected):
		assert expected == sub_gvcf._get_frame_for_genomic_position(genome_position)

	@pytest.mark.parametrize("protein_name,codon_number,expected",[(None,None,None),('nsp1',1,'ATG'),('S',607,'CAG'),('ORF10',39,'TAG')])
	def test_get_protein_codon_DNA_seq(self,sub_gvcf,protein_name,codon_number,expected):
		assert expected == sub_gvcf._get_protein_codon_DNA_seq(protein_name,codon_number)

	@pytest.mark.parametrize("genome_position1,genome_position2,expected",[(
		100,200,None),(300,900,None),(269,271,'E'),(270,271,'E'),(22369,22376,'YLQP'),(29558,29674,'MGYINVFAFPFTIYSLLLCRMNSRNYIAQVDVVNFNLT*')
	])
	def test_get_ref_AA_full_codon_seqs_between_two_genome_positions(self,sub_gvcf,genome_position1,genome_position2,expected):
		assert expected == sub_gvcf._get_ref_AA_full_codon_seqs_between_two_genome_positions(genome_position1,genome_position2)
