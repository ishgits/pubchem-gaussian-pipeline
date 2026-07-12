%nprocshared=16
%chk=water_F.chk
# opt=(tight,calcfc) b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water)

water PCM 298 K 6-311++G(2df,2p)

0 1
O        1.08126000   0.03580000  -0.02748000
H        0.80241000   0.77769000  -0.58499000
H        2.04914000   0.07292000  -0.05537000

--Link1--
%nprocshared=16
%chk=water_F.chk
# freq b3lyp/6-311++g(2df,2p) scrf=(iefpcm,solvent=water) temperature=298 Geom=AllChk Guess=Read
