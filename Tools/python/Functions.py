#!/usr/bin/env python
'''@package docstring
Loads some numerical functions defined in a C++ file
'''

from PandaCore.Utils.load import Load 
from PandaCore.Utils.root import root 

Load('Functions')

# need to instantiate these things
root.bound
root.dsign
root.Mxx
root.MT
root.SignedDeltaPhi
root.DeltaPhi
root.DeltaR2
root.ExpErf
root.metsf
root.metsf2018
root.fr_pho
root.sfpho
root.sfphonew
root.sfpho2018
root.photrigsf
root.photrigsf2018
root.combinedPhi
root.combinedPt
root.Expnom
root.sfwjets
root.sfgjets
root.sfzjets
root.sfzvv
root.tightIDsf
root.eleunc
root.tauSF
