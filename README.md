# QED package

**QED** stands for quantitative estimation of drug-likeness and the concept has been introduced by Richard Bickerton and coworkers [Bickerton, G.R.; Paolini, G.V.; Besnard, J.; Muresan, S.; Hopkins, A.L. (2012) ‘Quantifying the chemical beauty of drugs’, Nature Chemistry, 4, 90-98](http://dx.doi.org/10.1038/nchem.1243). This module relies on [RDKit](https://www.rdkit.org) as a chem-informatics toolkit.

## Introduction

This section is about installing and using **QED** within your own Python scripts or as a standalone Python tool.

The empirical rationale of the **QED** measure reflects the underlying distribution of molecular properties including molecular weight, logP, topological polar surface area, number of hydrogen bond donors and acceptors, the number of aromatic rings and rotatable bonds, and the presence of unwanted chemical functionalities.

The **QED** results as generated by this module are not completely identical to those from the [original publication](http://dx.doi.org/10.1038/nchem.1243). These differences are a consequence of differences within the underlying calculated property calculators used in both methods. For example, discrepancies can be noted in the results from the logP calculations, nevertheless despite the fact that both approaches (Pipeline Pilot in the original publication and RDKit in this implementation) mention to use the [Wildmann and Crippen](http://pubs.acs.org/doi/abs/10.1021/ci990307l) methodology for the calculation of their logP-values. In this respect, Gregory Gerebtzoff has been so kind to perform a refitting of the **QED** parameters with logP values generated by RDKit ([see rdkit-discuss](http://sourceforge.net/mailarchive/forum.php?thread_name=E05E80C886E33E4BA10E5686C606617602FF1D7A14%40RKAMSEM707.emea.roche.com&forum_name=rdkit-discuss)). These refitted values have been implemented in **QED** as the default values; however, the original publication values can still be used if desired.
This section assumes you have installed RDKit correctly and that you are familiar with the basic functions of it. It is also recommended to have read the original [QED publication](http://dx.doi.org/10.1038/nchem.1243).

## Download and installation

Download **QED** from GitHub (in this section we assume you have downloaded the file into your `~/Downloads` directory) and un-tar this file into this directory:

```console
> cd ~/Downloads
> sudo tar -xvf qed-1.0.1.tar.gz
> cd qed-1.0.1
```

You should now have a number of files in your `~/Downloads/qed-1.0.1` directory:

```console
> ls -l
qed-1.0.1/
qed-1.0.1/LICENSE
qed-1.0.1/README.md
qed-1.0.1/dist/
qed-1.0.1/dist/qed-1.0.1.tar.gz
qed-1.0.1/how_to_make_distribution.txt
qed-1.0.1/qed/
qed-1.0.1/qed/__init__.py
qed-1.0.1/qed/qed.py
qed-1.0.1/setup.py
```

Move into the `dist` directory and `untar` again:

```
> cd dist
> tar -xvf qed-1.0.1.tar.gz
> cd qed-1.0.1
```

Now install with Python:

```console
> python setup.py install
```

This process creates a `qed` folder with all the `qed` files into your default Python site-package install directory. It might be that you need to get root permissions:

```console
> sudo python setup.py install
```

Check your installation by launching `Python`:

```python
>>> from qed import qed
>>> from rdkit import Chem
>>> m = Chem.MolFromSmiles('c1ccccc1')
>>> qed.default(m)
0.44619898311523704
```

## Using QED

The `qed()` function takes as argument a RDKit molecule and returns the corresponding QED value calculated from it.
The `qed()` function comes in three flavors, each differing in the relative weight that is imposed on the underlying [molecular descriptors](http://dx.doi.org/10.1038/nchem.1243). These three flavors correspond to the three different **QED** measures that were described in the original publication:

- QED<sub>w,max</sub> using the set of weights that give maximal information content.
- QED<sub>w,mo</sub> using the mean weights of the optimal 1,000 weight combinations that give the highest information content.
- QED<sub>w,u</sub> with all weights as unity, hence unweighted.

Specifying the required QED weighting scheme in **QED** is done using the corresponding function:

- `qed.weights_mean()` uses the mean weighting scheme and corresponds to QED<sub>w,mo</sub>. Another name for this function is `qed.default()`.
- `qed.weights_max()` specifies the max weighting scheme and corresponds to QED<sub>w,max</sub>.
- `qed.weights_none()` specifies unit weights and corresponds to QED<sub>w,u</sub>.

and exemplified below:

```python
>>> from qed import qed
>>> from rdkit import Chem
>>> m = Chem.MolFromSmiles('c1ccccc1')
>>> qed.default(m)
0.44619898311523704
>>> qed.weights_mean(m)
0.44619898311523704
>>> qed.weights_max(m)
0.4733526950948539
>>> qed.weights_none(m)
0.3047153431243204
```

As already mentioned above, the current implementation of **QED** uses the refitted logP parameters from Gregory Gerebtzoff. However, the original values can still be used by specifying `False` as second argument to the appropriate function call:

```python
>>> qed.default(m, False)
0.4426283718993647
>>> qed.weights_mean(m, False)
0.4426283718993647
>>> qed.weights_max(m, False)
0.4706596045646091
>>> qed.weights_none(m, False)
0.3021185764176498
```

The **QED** module also provides its own test set:

```python
>>> test = qed.test_data()
>>> for name in test: print(test[name], name)
...
Nc1nc(NC2CC2)c2ncn([C@@H]3C[C@H](CO)C=C3)c2n1 Abacavir
CC(=O)NCCCS(O)(=O)=O Acamprosate
CCCC(=O)Nc1ccc(OCC(O)CNC(C)C)c(c1)C(C)=O Acebutolol
CC(=O)Nc1ccc(O)cc1 Acetaminophen
CC(=O)Nc1nnc(s1)S(N)(=O)=O Acetazolamide
CC(=O)c1ccc(cc1)S(=O)(=O)NC(=O)NC1CCCCC1 Acetohexamide
CC(=O)c1ccc2Sc3ccccc3N(CCCN3CCN(CCO)CC3)c2c1 Acetophenazine
Fc4ccc(C1CCNCC1COc3ccc2OCOc2c3)cc4 Paroxetine
Cc1oncc1C(=O)Nc2ccc(C(F)(F)F)cc2 Leflunomide
CN1C4CCCC1CC(NC(=O)c2nn(C)c3ccccc23)C4 Granisetron
CCCN2CC(CSC)CC1c3cccc4[nH]cc(CC12)c34 Pergolide
CCc3c(C)[nH]c2CCC(CN1CCOCC1)C(=O)c23 Molindone
CCCCCCCCCCCCCCCC(=O)OCC(NC(=O)C(Cl)Cl)C(O)c1ccc([N+]([O-])=O)cc1 ChloramphenicalPalmitate
CCCCCCCCCCCCCCCOC(=O)C2C(O)C(O)C(C(NC(=O)C1CC(CCC)CN1C)C(C)Cl)OC2SC ClindamycinPalmitate
CCOc3nc2cccc(C(=O)OC(C)OC(=O)OC1CCCCC1)c2n3Cc6ccc(c4ccccc4c5nn[nH]n5)cc6 CandesartanCilexetil
CN(C)CCC=c2c1ccccc1sc3ccc(Cl)cc23 Chlorprothixene
O=c3c(O)c(C2CCC(c1ccc(Cl)cc1)CC2)c(=O)c4ccccc34 Atovaquone
CN(C)CCCN3c1ccccc1CCc2ccc(Cl)cc23 Clomipramine
CN4CCCC(CC3c1ccccc1Sc2ccccc23)C4 Methixene
CCN(CC)C(C)Cn3c1ccccc1sc2ccccc23 Ethopropazine
N=C(CCSCc1csc(N=C(N)N)n1)NS(N)(=O)=O Famotidine
CNC(=NCCSCc1nc[nH]c1C)NC#N Cimetidine
CCCCCNC(=N)NN=Cc1c[nH]c2ccc(CO)cc12 Tegaserod
C=CC3=C(C(=O)O)N2C(=O)C(NC(=O)C(=NO)c1csc(N)n1)C2SC3 Cefdinir
CC5(C)SC4C(NC(=O)C(C(=O)Oc2ccc1CCCc1c2)c3ccccc3)C(=O)N4C5C(=O)O CarbenicillinIndanyl
```

One can inspect the individual properties that are used to calculate the Qed value by calling the `properties()` function:

```python
>>> for name in test:
...     mol = Chem.MolFromSmiles(test[name])
...     p = qed.properties(mol)
...     print("%6.2f\t%6.3f\t%6d\t%6d\t%6.2f\t%6d\t%6d\t%6d\t%6.3f\t%-s" % (p[0],p[1],p[2],p[3],p[4],p[5],p[6],p[7],qed.default(mol),name))
...
286.34	 1.092	     6	     3	101.88	     4	     2	     1	 0.737	Abacavir
181.21	-0.600	     4	     2	 83.47	     4	     0	     2	 0.467	Acamprosate
336.43	 2.365	     5	     3	 87.66	    10	     1	     1	 0.571	Acebutolol
151.16	 1.351	     2	     2	 49.33	     1	     1	     1	 0.602	Acetaminophen
222.25	-0.856	     5	     2	115.04	     2	     1	     1	 0.662	Acetazolamide
324.40	 2.210	     4	     2	 92.34	     4	     1	     1	 0.833	Acetohexamide
411.57	 3.492	     6	     1	 47.02	     7	     2	     1	 0.688	Acetophenazine
329.37	 3.327	     4	     1	 39.72	     4	     2	     0	 0.917	Paroxetine
270.21	 3.254	     3	     1	 55.13	     2	     2	     0	 0.896	Leflunomide
312.42	 2.318	     3	     1	 50.16	     2	     2	     0	 0.927	Granisetron
314.50	 4.271	     2	     1	 19.03	     4	     2	     0	 0.871	Pergolide
276.38	 1.963	     3	     1	 45.33	     3	     1	     0	 0.923	Molindone
561.55	 6.941	     6	     2	118.77	    21	     1	     5	 0.056	ChloramphenicalPalmitate
663.41	 6.279	     8	     3	108.33	    22	     0	     3	 0.071	ClindamycinPalmitate
610.67	 6.319	    10	     1	143.34	    10	     5	     2	 0.141	CandesartanCilexetil
315.87	 5.188	     2	     0	  3.24	     3	     3	     0	 0.629	Chlorprothixene
366.84	 5.505	     3	     1	 54.37	     2	     2	     0	 0.741	Atovaquone
314.86	 4.528	     2	     0	  6.48	     4	     2	     0	 0.782	Clomipramine
309.48	 5.015	     2	     0	  3.24	     2	     2	     0	 0.735	Methixene
312.48	 5.020	     3	     0	  6.48	     5	     2	     0	 0.734	Ethopropazine
337.46	-0.558	     6	     5	173.33	     7	     1	     3	 0.263	Famotidine
252.35	 0.597	     5	     3	 88.89	     5	     1	     5	 0.239	Cimetidine
301.39	 2.298	     3	     5	 96.29	     7	     2	     4	 0.235	Tegaserod
395.42	-0.172	     8	     4	158.21	     5	     1	     4	 0.239	Cefdinir
494.57	 2.496	     7	     2	113.01	     6	     2	     4	 0.274	CarbenicillinIndanyl
>>> print("    MW\t ALOGP\t   HBA\t   HBD\t   PSA\t  ROTB\t  AROM\tALERTS\t   QED\tNAME")
    MW	 ALOGP	   HBA	   HBD	   PSA	  ROTB	  AROM	ALERTS	   QED	NAME
```

## Revision history

### Version 1.0.1
Incorporation of the refitted logP parameters of Gregory Gerebtzoff and making these values default ([rdkit-discuss](http://sourceforge.net/mailarchive/forum.php?thread_name=E05E80C886E33E4BA10E5686C606617602FF1D7A14%40RKAMSEM707.emea.roche.com&forum_name=rdkit-discuss)).
Modification of the `qed()` function to enable the selection of the original parameters, if desired.
