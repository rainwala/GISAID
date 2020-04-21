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

	@pytest.fixture
	def ins_gvcf2(self):
		ins_gvcf_line = GVCFLine.from_line("NC_045512\t.\t22304\tT\tTCCCACCAGA\n")
		return GisVar(ins_gvcf_line,['EPI_ISL_420726'])

	@pytest.fixture
	def ins_gvcf3(self):
		ins_gvcf_line = GVCFLine.from_line("NC_045512\t.\t54\tT\tTC\n")
		return GisVar(ins_gvcf_line,['EPI_ISL_416654'])

	@pytest.fixture
	def ins_gvcf4(self):
		ins_gvcf_line = GVCFLine.from_line("NC_045512\t.\t21569\tG\tGT\n")
		return GisVar(ins_gvcf_line,['EPI_ISL_420942'])

	@pytest.fixture
	def sub_gvcf_list(self,sub_gvcf,sub_gvcf2,sub_gvcf3):		
		return [sub_gvcf,sub_gvcf2,sub_gvcf3]

	@pytest.fixture
	def del_gvcf_list(self,del_gvcf,del_gvcf2,del_gvcf3,del_gvcf4,del_gvcf5):
		return [del_gvcf,del_gvcf2,del_gvcf3,del_gvcf4,del_gvcf5]

	@pytest.fixture
	def ins_gvcf_list(self,ins_gvcf,ins_gvcf2,ins_gvcf3,ins_gvcf4):
		return [ins_gvcf,ins_gvcf2,ins_gvcf3,ins_gvcf4]

	@pytest.fixture
	def all_gvcf_dict(self,sub_gvcf_list,del_gvcf_list,ins_gvcf_list):
		return {
			'sub': sub_gvcf_list,
			'del': del_gvcf_list,
			'ins': ins_gvcf_list
		}		

	""" GisVar level tests """

	def test_get_type(self,sub_gvcf,del_gvcf,ins_gvcf):
		assert 'sub' == sub_gvcf.get_type()
		assert 'del' == del_gvcf.get_type()
		assert 'ins' == ins_gvcf.get_type()

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

	@pytest.mark.parametrize("index,expected", [(0,'NC_045512.2:g.21590C>T'),(1,'NC_045512.2:g.241C>T'),(2,'NC_045512.2:g.1059C>T')])
	def test_sub_gvcf_make_genomic_variant_name(self,sub_gvcf_list,index,expected):
		assert expected == sub_gvcf_list[index].get_var_type_child().make_genomic_variant_name()

	@pytest.mark.parametrize("index,expected", [(0,'S:p.Leu10Leu'),(1,None),(2,'nsp2:p.Thr85Ile')])
	def test_sub_gvcf_make_protein_variant_name(self,sub_gvcf_list,index,expected):
		assert expected == sub_gvcf_list[index].get_var_type_child().make_protein_variant_name()

	@pytest.mark.parametrize("index,expected", [
		(0,'NC_045512.2:g.21990_21992del'),(1,'NC_045512.2:g.21991_21993del'),(2,'NC_045512.2:g.263del'),
		(3,'NC_045512.2:g.27426_27443del'), (4,'NC_045512.2:g.22928del')
	])
	def test_del_gvcf_make_genomic_variant_name(self,del_gvcf_list,index,expected):
		assert expected == del_gvcf_list[index].get_var_type_child().make_genomic_variant_name()
	
	@pytest.mark.parametrize("index,expected", [(0,'S:p.Val143_Tyr144delinsAsp'),(1,'S:p.Tyr144del'),(2,None),(3,'ORF7a:p.Leu12_Leu17del'),(4,'S:p.(Phe456fs)')])
	def test_del_gvcf_make_protein_variant_name(self,del_gvcf_list,index,expected):
		assert expected == del_gvcf_list[index].get_var_type_child().make_protein_variant_name()

	@pytest.mark.parametrize("index,expected", [
		(0,'NC_045512.2:g.3336_3337insATC'),(1,'NC_045512.2:g.22304_22305insCCCACCAGA'),
		(2,'NC_045512.2:g.54_55insC'),(3,'NC_045512.2:g.21569_21570insT')
	])
	def test_ins_gvcf_make_genomic_variant_name(self,ins_gvcf_list,index,expected):
		assert expected == ins_gvcf_list[index].get_var_type_child().make_genomic_variant_name()

	@pytest.mark.parametrize("index,expected", [(0,'nsp3:p.Glu206_Val207insSer'),(1,'S:p.Tyr248delinsSerHisGlnAsn'),(2,None),(3,'S:p.(Val3fs)')])
	def test_ins_gvcf_make_protein_variant_name(self,ins_gvcf_list,index,expected):
		assert expected == ins_gvcf_list[index].get_var_type_child().make_protein_variant_name()

	@pytest.mark.parametrize("type,index,expected", [
		('sub',0,True),('sub',1,True),('sub',2,False),
		('del',0,False),('del',1,False),('del',2,True),('del',3,False),('del',4,False),
		('ins',0,False),('ins',1,False),('ins',2,True),('ins',3,False),
	])
	def test_is_synonymous(self,all_gvcf_dict,type,index,expected):
		assert expected == all_gvcf_dict[type][index].is_synonymous( all_gvcf_dict[type][index].get_var_type_child().make_protein_variant_name() )
