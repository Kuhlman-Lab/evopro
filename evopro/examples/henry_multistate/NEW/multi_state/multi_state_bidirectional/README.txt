# score function is simple AF2 confidence (pLDDT) but averaged for both conformations
# AF2 predictions are done SEPARATELY but MPNN predictions are TOGETHER

# here we reverse the sequence symmetry using the preprocessing script generate_json_multi_state.py
# then we apply bidirectional coding filter using ProteinMPNN with the --bidirectional flag in evopro
