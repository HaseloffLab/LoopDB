from loopDB import LoopDB

loopDB = LoopDB( "postgresql:///loopdb", clean = True )

BsaI = loopDB.addRE(name = "BsaI",	seq = "GGTCTCA")
SapI = loopDB.addRE(name = "SapI",	seq = "GGTCTCA")

AB = loopDB.addRES(name = "AB",	re = BsaI,	site5 = "AAAA",	site3 = "BBBB")
BC = loopDB.addRES(name = "BC",	re = BsaI,	site5 = "BBBB",	site3 = "CCCC")
AC = loopDB.addRES(name = "AC",	re = BsaI,	site5 = "AAAA",	site3 = "CCCC")

XY = loopDB.addRES(name = "XY",	re = SapI,	site5 = "XXXX",	site3 = "YYYY")
YZ = loopDB.addRES(name = "YZ",	re = SapI,	site5 = "YYYY",	site3 = "ZZZZ")
XZ = loopDB.addRES(name = "XZ",	re = SapI,	site5 = "XXXX",	site3 = "ZZZZ")

L0 = loopDB.addBaseSeq(name = "L0",	seq = "CCCCCCCCCCCCCCC")
Ly = loopDB.addBaseSeq(name = "LY",	seq = "GGGGGGGGGGGGGGG",	receiver = AC)
Lx = loopDB.addBaseSeq(name = "LX",	seq = "TTTTTTTTTTTTTTT",	receiver = XZ)

Ly1 = loopDB.addBackbone(name = "Ly1",	baseSeq = Ly,	adapter = XY)
Ly2 = loopDB.addBackbone(name = "Ly2",	baseSeq = Ly,	adapter = YZ)                                        

Lx1 = loopDB.addBackbone(name = "Lx1",	baseSeq = Lx,	adapter = AB)
Lx2 = loopDB.addBackbone(name = "Lx2",	baseSeq = Lx,	adapter = BC)

L0A  = loopDB.addBackbone(name = "L0A",	baseSeq = L0,	adapter = AB)
L0B  = loopDB.addBackbone(name = "L0B",	baseSeq = L0,	adapter = BC)

part01 = loopDB.addPart(name = "Part01",  seq = "ATATAT", backbone = L0A)
part02 = loopDB.addPart(name = "Part02",  seq = "GCGCGC", backbone = L0B)

part03 = loopDB.addPart(name = "Part03",  seq = "TCTCTC", backbone = L0A)
part04 = loopDB.addPart(name = "Part04",  seq = "AGAGAG", backbone = L0B)

part11 = loopDB.addPart(name = "Part11", backbone = Ly1, children = [part01, part02])
part12 = loopDB.addPart(name = "Part12", backbone = Ly2, children = [part03, part04])

feature = loopDB.addFeature(start = 10, end = 20, forward = True, label ="Feature", type = "misc")
part20 = loopDB.addPart(name = "Part20",   children = [part11, part12], backbone = Lx1, features = [feature])

print part20.fullSeq

loopDB.commit()