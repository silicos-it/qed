# Copyright 2012 by Silicos-it, a division of Imacosi BVBA
# Copyright 2020- by UAMC, the Medicinal Chemistry department of the University
# of Antwerp


__all__ = ['weights_max', 'weights_mean', 'weights_none', 'default']
__version__ = "1.0.1"


# RDKit
from rdkit.Chem import Descriptors
from rdkit import Chem

# General
from copy import deepcopy
from math import exp, log
import sys, os


#
AliphaticRings = Chem.MolFromSmarts('[$([A;R][!a])]')


#
AcceptorSmarts = [
	'[oH0;X2]',
	'[OH1;X2;v2]',
	'[OH0;X2;v2]',
	'[OH0;X1;v2]',
	'[O-;X1]',
	'[SH0;X2;v2]',
	'[SH0;X1;v2]',
	'[S-;X1]',
	'[nH0;X2]',
	'[NH0;X1;v3]',
	'[$([N;+0;X3;v3]);!$(N[C,S]=O)]'
	]
Acceptors = []
for hba in AcceptorSmarts:
	Acceptors.append(Chem.MolFromSmarts(hba))



#
StructuralAlertSmarts = [
	'*1[O,S,N]*1',
	'[S,C](=[O,S])[F,Br,Cl,I]',
	'[CX4][Cl,Br,I]',
	'[C,c]S(=O)(=O)O[C,c]',
	'[$([CH]),$(CC)]#CC(=O)[C,c]',
	'[$([CH]),$(CC)]#CC(=O)O[C,c]',
	'n[OH]',
	'[$([CH]),$(CC)]#CS(=O)(=O)[C,c]',
	'C=C(C=O)C=O',
	'n1c([F,Cl,Br,I])cccc1',
	'[CH1](=O)',
	'[O,o][O,o]',
	'[C;!R]=[N;!R]',
	'[N!R]=[N!R]',
	'[#6](=O)[#6](=O)',
	'[S,s][S,s]',
	'[N,n][NH2]',
	'C(=O)N[NH2]',
	'[C,c]=S',
	'[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]=[$([CH2]),$([CH][CX4]),$(C([CX4])[CX4])]',
	'C1(=[O,N])C=CC(=[O,N])C=C1',
	'C1(=[O,N])C(=[O,N])C=CC=C1',
	'a21aa3a(aa1aaaa2)aaaa3',
	'a31a(a2a(aa1)aaaa2)aaaa3',
	'a1aa2a3a(a1)A=AA=A3=AA=A2',
	'c1cc([NH2])ccc1',
	'[Hg,Fe,As,Sb,Zn,Se,se,Te,B,Si,Na,Ca,Ge,Ag,Mg,K,Ba,Sr,Be,Ti,Mo,Mn,Ru,Pd,Ni,Cu,Au,Cd,Al,Ga,Sn,Rh,Tl,Bi,Nb,Li,Pb,Hf,Ho]',
	'I',
	'OS(=O)(=O)[O-]',
	'[N+](=O)[O-]',
	'C(=O)N[OH]',
	'C1NC(=O)NC(=O)1',
	'[SH]',
	'[S-]',
	'c1ccc([Cl,Br,I,F])c([Cl,Br,I,F])c1[Cl,Br,I,F]',
	'c1cc([Cl,Br,I,F])cc([Cl,Br,I,F])c1[Cl,Br,I,F]',
	'[CR1]1[CR1][CR1][CR1][CR1][CR1][CR1]1',
	'[CR1]1[CR1][CR1]cc[CR1][CR1]1',
	'[CR2]1[CR2][CR2][CR2][CR2][CR2][CR2][CR2]1',
	'[CR2]1[CR2][CR2]cc[CR2][CR2][CR2]1',
	'[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1',
	'[CH2R2]1N[CH2R2][CH2R2][CH2R2][CH2R2][CH2R2][CH2R2]1',
	'C#C',
	'[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]@[CR2]@[CR2]@[OR2,NR2]',
	'[$([N+R]),$([n+R]),$([N+]=C)][O-]',
	'[C,c]=N[OH]',
	'[C,c]=NOC=O',
	'[C,c](=O)[CX4,CR0X3,O][C,c](=O)',
	'c1ccc2c(c1)ccc(=O)o2',
	'[O+,o+,S+,s+]',
	'N=C=O',
	'[NX3,NX4][F,Cl,Br,I]',
	'c1ccccc1OC(=O)[#6]',
	'[CR0]=[CR0][CR0]=[CR0]',
	'[C+,c+,C-,c-]',
	'N=[N+]=[N-]',
	'C12C(NC(N1)=O)CSC2',
	'c1c([OH])c([OH,NH2,NH])ccc1',
	'P',
	'[N,O,S]C#N',
	'C=C=O',
	'[Si][F,Cl,Br,I]',
	'[SX2]O',
	'[SiR0,CR0](c1ccccc1)(c2ccccc2)(c3ccccc3)',
	'O1CCCCC1OC2CCC3CCCCC3C2',
	'N=[CR0][N,n,O,S]',
	'[cR2]1[cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2][cR2]1[cR2]2[cR2][cR2][cR2]([Nv3X3,Nv4X4])[cR2][cR2]2',
	'C=[C!r]C#N',
	'[cR2]1[cR2]c([N+0X3R0,nX3R0])c([N+0X3R0,nX3R0])[cR2][cR2]1',
	'[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2]c([N+0X3R0,nX3R0])[cR2]1',
	'[cR2]1[cR2]c([N+0X3R0,nX3R0])[cR2][cR2]c1([N+0X3R0,nX3R0])',
	'[OH]c1ccc([OH,NH2,NH])cc1',
	'c1ccccc1OC(=O)O',
	'[SX2H0][N]',
	'c12ccccc1(SC(S)=N2)',
	'c12ccccc1(SC(=S)N2)',
	'c1nnnn1C=O',
	's1c(S)nnc1NC=O',
	'S1C=CSC1=S',
	'C(=O)Onnn',
	'OS(=O)(=O)C(F)(F)F',
	'N#CC[OH]',
	'N#CC(=O)',
	'S(=O)(=O)C#N',
	'N[CH2]C#N',
	'C1(=O)NCC1',
	'S(=O)(=O)[O-,OH]',
	'NC[F,Cl,Br,I]',
	'C=[C!r]O',
	'[NX2+0]=[O+0]',
	'[OR0,NR0][OR0,NR0]',
	'C(=O)O[C,H1].C(=O)O[C,H1].C(=O)O[C,H1]',
	'[CX2R0][NX3R0]',
	'c1ccccc1[C;!R]=[C;!R]c2ccccc2',
	'[NX3R0,NX4R0,OR0,SX2R0][CX4][NX3R0,NX4R0,OR0,SX2R0]',
	'[s,S,c,C,n,N,o,O]~[n+,N+](~[s,S,c,C,n,N,o,O])(~[s,S,c,C,n,N,o,O])~[s,S,c,C,n,N,o,O]',
	'[s,S,c,C,n,N,o,O]~[nX3+,NX3+](~[s,S,c,C,n,N])~[s,S,c,C,n,N]',
	'[*]=[N+]=[*]',
	'[SX3](=O)[O-,OH]',
	'N#N',
	'F.F.F.F',
	'[R0;D2][R0;D2][R0;D2][R0;D2]',
	'[cR,CR]~C(=O)NC(=O)~[cR,CR]',
	'C=!@CC=[O,S]',
	'[#6,#8,#16][C,c](=O)O[C,c]',
	'c[C;R0](=[O,S])[C,c]',
	'c[SX2][C;!R]',
	'C=C=C',
	'c1nc([F,Cl,Br,I,S])ncc1',
	'c1ncnc([F,Cl,Br,I,S])c1',
	'c1nc(c2c(n1)nc(n2)[F,Cl,Br,I])',
	'[C,c]S(=O)(=O)c1ccc(cc1)F',
	'[15N]',
	'[13C]',
	'[18O]',
	'[34S]'
	]
StructuralAlerts = []
for smarts in StructuralAlertSmarts:
	StructuralAlerts.append(Chem.MolFromSmarts(smarts))



# ADS parameters for the 8 molecular properties: [row][column]
# 	rows[8]: 	MW, ALOGP, HBA, HBD, PSA, ROTB, AROM, ALERTS
# 	columns[7]: A, B, C, D, E, F, DMAX
# ALOGP parameters from Gregory Gerebtzoff (2012, Roche)
pads1 = [	[2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561],
			[0.486849448, 186.2293718, 2.066177165, 3.902720615, 1.027025453, 0.913012565, 145.4314800],
			[2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046],
			[1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616],
			[1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167],
			[0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403],
			[3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610],
			[0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140]		]
# ALOGP parameters from the original publication
pads2 = [	[2.817065973, 392.5754953, 290.7489764, 2.419764353, 49.22325677, 65.37051707, 104.9805561],
			[3.172690585, 137.8624751, 2.534937431, 4.581497897, 0.822739154, 0.576295591, 131.3186604],
			[2.948620388, 160.4605972, 3.615294657, 4.435986202, 0.290141953, 1.300669958, 148.7763046],
			[1.618662227, 1010.051101, 0.985094388, 0.000000001, 0.713820843, 0.920922555, 258.1632616],
			[1.876861559, 125.2232657, 62.90773554, 87.83366614, 12.01999824, 28.51324732, 104.5686167],
			[0.010000000, 272.4121427, 2.558379970, 1.565547684, 1.271567166, 2.758063707, 105.4420403],
			[3.217788970, 957.7374108, 2.274627939, 0.000000001, 1.317690384, 0.375760881, 312.3372610],
			[0.010000000, 1199.094025, -0.09002883, 0.000000001, 0.185904477, 0.875193782, 417.7253140]		]



def ads(x, a, b, c, d, e, f, dmax):
	return ((a+(b/(1+exp(-1*(x-c+d/2)/e))*(1-1/(1+exp(-1*(x-c-d/2)/f))))) / dmax)



def properties(mol):
	"""
	Calculates the properties that are required to calculate the QED descriptor.
	"""
	matches = []
	if (mol is None):
		print("Error: in function properties(mol): mol argument is \'None\'")
		sys.exit(1)
	x = [0] * 8
	x[0] = Descriptors.MolWt(mol)												# MW
	x[1] = Descriptors.MolLogP(mol)												# ALOGP
	for hba in Acceptors:														# HBA
		if (mol.HasSubstructMatch(hba)):
			matches = mol.GetSubstructMatches(hba)
			x[2] += len(matches)
	x[3] = Descriptors.NumHDonors(mol)											# HBD
	x[4] = Descriptors.TPSA(mol)												# PSA
	x[5] = Descriptors.NumRotatableBonds(mol)									# ROTB
	x[6] = Chem.GetSSSR(Chem.DeleteSubstructs(deepcopy(mol), AliphaticRings))	# AROM
	for alert in StructuralAlerts:												# ALERTS
		if (mol.HasSubstructMatch(alert)): x[7] += 1
	return x



def qed(w, p, gerebtzoff):
	d = [0.00] * 8
	if (gerebtzoff):
		for i in range(0, 8):
			d[i] = ads(p[i], pads1[i][0], pads1[i][1], pads1[i][2], pads1[i][3], pads1[i][4], pads1[i][5], pads1[i][6])
	else:
		for i in range(0, 8):
			d[i] = ads(p[i], pads2[i][0], pads2[i][1], pads2[i][2], pads2[i][3], pads2[i][4], pads2[i][5], pads2[i][6])
	t = 0.0
	for i in range(0, 8):
		t += w[i] * log(d[i])
	return (exp(t / sum(w)))



def weights_max(mol, gerebtzoff = True):
	"""
	Calculates the QED descriptor using maximal descriptor weights.
	"""
	props = properties(mol)
	return qed([0.50, 0.25, 0.00, 0.50, 0.00, 0.50, 0.25, 1.00], props, gerebtzoff)



def weights_mean(mol, gerebtzoff = True):
	"""
	Calculates the QED descriptor using average descriptor weights.
	"""
	props = properties(mol)
	return qed([0.66, 0.46, 0.05, 0.61, 0.06, 0.65, 0.48, 0.95], props, gerebtzoff)



def weights_none(mol, gerebtzoff = True):
	"""
	Calculates the QED descriptor using unit weights.
	"""
	props = properties(mol)
	return qed([1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00], props, gerebtzoff)



def default(mol, gerebtzoff = True):
	"""
	Calculates the QED descriptor using average descriptor weights and Gregory Gerebtzoff parameters.
	"""
	return weights_mean(mol, gerebtzoff)



def test_data():
	"""
	Returns some test data
	"""
	return({
	"Abacavir" : "Nc1nc(NC2CC2)c2ncn([C@@H]3C[C@H](CO)C=C3)c2n1",
	"Acamprosate" : "CC(=O)NCCCS(O)(=O)=O",
	"Acebutolol" : "CCCC(=O)Nc1ccc(OCC(O)CNC(C)C)c(c1)C(C)=O",
	"Acetaminophen" : "CC(=O)Nc1ccc(O)cc1",
	"Acetazolamide" : "CC(=O)Nc1nnc(s1)S(N)(=O)=O",
	"Acetohexamide" : "CC(=O)c1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1",
	"Acetophenazine" : "CC(=O)c1ccc2Sc3ccccc3N(CCCN3CCN(CCO)CC3)c2c1",
	"Paroxetine" : "Fc4ccc(C1CCNCC1COc3ccc2OCOc2c3)cc4",
	"Leflunomide" : "Cc1oncc1C(=O)Nc2ccc(C(F)(F)F)cc2",
	"Granisetron" : "CN1C4CCCC1CC(NC(=O)c2nn(C)c3ccccc23)C4",
	"Pergolide" : "CCCN2CC(CSC)CC1c3cccc4[nH]cc(CC12)c34",
	"Molindone" : "CCc3c(C)[nH]c2CCC(CN1CCOCC1)C(=O)c23",
	"ChloramphenicalPalmitate" : "CCCCCCCCCCCCCCCC(=O)OCC(NC(=O)C(Cl)Cl)C(O)c1ccc([N+]([O-])=O)cc1",
	"ClindamycinPalmitate" : "CCCCCCCCCCCCCCCOC(=O)C2C(O)C(O)C(C(NC(=O)C1CC(CCC)CN1C)C(C)Cl)OC2SC",
	"CandesartanCilexetil" : "CCOc3nc2cccc(C(=O)OC(C)OC(=O)OC1CCCCC1)c2n3Cc6ccc(c4ccccc4c5nn[nH]n5)cc6",
 	"Chlorprothixene" : "CN(C)CCC=c2c1ccccc1sc3ccc(Cl)cc23",
	"Atovaquone" : "O=c3c(O)c(C2CCC(c1ccc(Cl)cc1)CC2)c(=O)c4ccccc34",
	"Clomipramine" : "CN(C)CCCN3c1ccccc1CCc2ccc(Cl)cc23",
	"Methixene" : "CN4CCCC(CC3c1ccccc1Sc2ccccc23)C4",
	"Ethopropazine" : "CCN(CC)C(C)Cn3c1ccccc1sc2ccccc23",
	"Famotidine" : "N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O",
	"Cimetidine" : "CNC(=NCCSCc1nc[nH]c1C)NC#N",
	"Tegaserod" : "CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12",
	"Cefdinir" : "C=CC3=C(C(=O)O)N2C(=O)C(NC(=O)C(=NO)c1csc(N)n1)C2SC3",
	"CarbenicillinIndanyl" : "CC5(C)SC4C(NC(=O)C(C(=O)Oc2ccc1CCCc1c2)c3ccccc3)C(=O)N4C5C(=O)O",
	})
