from gisaid_variant import GisVar,GisSub,GisDel,GisIns
from gisaid_record import GVCFLine
import pytest

class TestGisVar:

	""" the following are useful fixtures to pass around for various tests """
	@pytest.fixture
	def sub_gvcf(self):
		sub_gvcf_line = GVCFLine.from_line("NC_045512\t.\t21590\tC\tT\n")
		return GisVar(sub_gvcf_line,['EPI_ISL_421184'])

	@pytest.fixture
	def sub_gvcf2(self):
		sub_gvcf_line = GVCFLine.from_line("NC_045512\t.\t241\tC\tT\n")
		return GisVar(sub_gvcf_line,['EPI_ISL_421184'])

	@pytest.fixture
	def sub_gvcf3(self):
		sub_gvcf_line = GVCFLine.from_line("NC_045512\t.\t1059\tC\tT\n")
		return GisVar(sub_gvcf_line,['EPI_ISL_421184'])

	@pytest.fixture
	def del_gvcf(self):
		del_gvcf_line = GVCFLine.from_line("NC_045512\t.\t21989\tGTTT\tG\n")
		return GisVar(del_gvcf_line,['EPI_ISL_421184'])

	@pytest.fixture
	def del_gvcf2(self):
		del_gvcf_line = GVCFLine.from_line("NC_045512\t.\t21990\tTTTA\tT\n")
		return GisVar(del_gvcf_line,['EPI_ISL_413522'])

	@pytest.fixture
	def del_gvcf3(self):
		del_gvcf_line = GVCFLine.from_line("NC_045512\t.\t262\tTA\tT\n")
		return GisVar(del_gvcf_line,['EPI_ISL_414367'])

	@pytest.fixture
	def del_gvcf4(self):
		del_gvcf_line = GVCFLine.from_line("NC_045512\t.\t27425\tCACTCGCTACTTGTGAGCT\tC\n")
		return GisVar(del_gvcf_line,['EPI_ISL_422281'])
	
	@pytest.fixture
	def del_gvcf5(self):
		del_gvcf_line = GVCFLine.from_line("NC_045512\t.\t22927\tGT\tG")
		return GisVar(del_gvcf_line,['EPI_ISL_416036'])

	@pytest.fixture
	def ins_gvcf(self):
		ins_gvcf_line = GVCFLine.from_line("NC_045512\t.\t3336\tA\tAATC\n")
		return GisVar(ins_gvcf_line,['EPI_ISL_410721'])

	""" GisVar level tests """

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
		100,200,None),(300,900,None),(269,271,'GAG'),(270,271,'GAG'),(22369,22376,'TATCTTCAACCT'),
		(29558,29674,'ATGGGCTATATAAACGTTTTCGCTTTTCCGTTTACGATATATAGTCTACTCTTGTGCAGAATGAATTCTCGTAACTACATAGCACAAGTAGATGTAGTTAACTTTAATCTCACATAG')
	])
	def test_get_ref_DNA_full_codon_seqs_between_two_genome_positions(self,sub_gvcf,genome_position1,genome_position2,expected):
		assert expected == sub_gvcf._get_ref_DNA_full_codon_seqs_between_two_genome_positions(genome_position1,genome_position2)

	""" GisSub, GisDel, GisIns level tests"""

	def test_sub_gvcf_make_genomic_variant_name(self,sub_gvcf,sub_gvcf2,sub_gvcf3):
		assert 'NC_045512.2:g.21590C>T' == sub_gvcf.get_var_type_child().make_genomic_variant_name()
		assert 'NC_045512.2:g.241C>T' == sub_gvcf2.get_var_type_child().make_genomic_variant_name()	
		assert 'NC_045512.2:g.1059C>T' == sub_gvcf3.get_var_type_child().make_genomic_variant_name()

	def test_sub_gvcf_make_protein_variant_name(self,sub_gvcf,sub_gvcf2,sub_gvcf3):
		assert 'S:p.Leu10Leu' == sub_gvcf.get_var_type_child().make_protein_variant_name()
		assert None == sub_gvcf2.get_var_type_child().make_protein_variant_name()
		assert 'nsp2:p.Thr85Ile' == sub_gvcf3.get_var_type_child().make_protein_variant_name()

	def test_del_gvcf_make_genomic_variant_name(self,del_gvcf,del_gvcf2,del_gvcf3,del_gvcf4,del_gvcf5):
		assert 'NC_045512.2:g.21990_21992del' == del_gvcf.get_var_type_child().make_genomic_variant_name()
		assert 'NC_045512.2:g.21991_21993del' == del_gvcf2.get_var_type_child().make_genomic_variant_name()
		assert 'NC_045512.2:g.263del' == del_gvcf3.get_var_type_child().make_genomic_variant_name()
		assert 'NC_045512.2:g.27426_27443del' == del_gvcf4.get_var_type_child().make_genomic_variant_name()
		assert 'NC_045512.2:g.22928del' == del_gvcf5.get_var_type_child().make_genomic_variant_name()
	
	def test_del_gvcf_make_protein_variant_name(self,del_gvcf,del_gvcf2,del_gvcf3,del_gvcf4,del_gvcf5):
		assert 'S:p.Val143_Tyr144delinsAsp' == del_gvcf.get_var_type_child().make_protein_variant_name()
		assert 'S:p.Tyr144del' == del_gvcf2.get_var_type_child().make_protein_variant_name()
		assert None == del_gvcf3.get_var_type_child().make_protein_variant_name()
		assert 'ORF7a:p.Leu12_Leu17del' == del_gvcf4.get_var_type_child().make_protein_variant_name()
		assert 'S:p.(Phe456fs)'
