%nprocshared=16
%chk=glycine_F.chk
# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

glycine PCM 298 K 6-311++G(2df,2p)

0 1
O        0.78940000  -1.25734000  -1.12870000
O        0.97286000  -1.12975000   1.09127000
N       -1.80090000  -0.37860000   1.23048000
C       -1.26841000  -1.12584000   0.08183000
C        0.25383000  -1.14539000   0.10476000
H       -1.61240000  -0.64572000  -0.83836000
H       -1.62237000  -2.16105000   0.10639000
H       -2.80352000  -0.53714000   1.30266000
H       -1.37895000  -0.74896000   2.08458000
H        1.75531000  -1.25254000  -0.96232000

--Link1--
%nprocshared=16
%chk=glycine_F.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read
