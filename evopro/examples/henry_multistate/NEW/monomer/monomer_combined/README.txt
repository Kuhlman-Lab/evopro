# same as monomer_RMSD_to_starting EXCEPT:
# 1. uses a combination of pLDDT and RMSD-to-starting for score function
# combined score function: total = (-pLDDT/10) + RMSD * 5 
# lower is better, as with all Evopro score functions