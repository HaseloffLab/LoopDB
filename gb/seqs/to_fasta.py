from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation, CompoundLocation
from Bio.Alphabet import IUPAC
from Bio import SeqIO
import sys
import os.path


for record in SeqIO.parse(sys.argv[1], "genbank"):
	print ">"+sys.argv[1].split(".")[0]
	print record.seq