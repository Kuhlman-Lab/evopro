ptms_dict = {"TPO":"T", "SEP":"S",  "MSE":"M", "PTR":"Y", "PCA":"E", "HY3":"P", "PHO":"X","ACE":"X", "NH2":"X", "UNK":"X"}
ptms = ["TPO", "SEP",  "MSE", "PTR", "PCA", "HY3", "PHO","ACE", "NH2", "UNK"]


num2aa_full=[
    'ALA','ARG','ASN','ASP','CYS',
    'GLN','GLU','GLY','HIS','ILE',
    'LEU','LYS','MET','PHE','PRO',
    'SER','THR','TRP','TYR','VAL',
    'UNK','MAS',
    'DA','DC','DG','DT','DX',
    'RA','RC','RG','RU','RX',
    'Al', 'As', 'Au', 'B',
    'Be', 'Br', 'C', 'Ca', 'Cl',
    'Co', 'Cr', 'Cu', 'F', 'Fe',
    'Hg', 'I', 'Ir', 'K', 'Li', 'Mg',
    'Mn', 'Mo', 'N', 'Ni', 'O',
    'Os', 'P', 'Pb', 'Pd', 'Pr',
    'Pt', 'Re', 'Rh', 'Ru', 'S',
    'Sb', 'Se', 'Si', 'Sn', 'Tb',
    'Te', 'U', 'W', 'V', 'Y', 'Zn',
    'ATM'
]

num2aa=[
    'A','R','N','D','C',
    'Q','E','G','H','I',
    'L','K','M','F','P',
    'S','T','W','Y','V',
    '-','X',
    'a','c','g','t','x',
    'b','d','h','u','y',
    'Al', 'As', 'Au', 'Bo',
    'Be', 'Br', 'Cb', 'Ca', 'Cl',
    'Co', 'Cr', 'Cu', 'Fl', 'Fe',
    'Hg', 'Io', 'Ir', 'Kt', 'Li', 'Mg',
    'Mn', 'Mo', 'Nk', 'Ni', 'Ox',
    'Os', 'Ph', 'Pb', 'Pd', 'Pr',
    'Pt', 'Re', 'Rh', 'Ru', 'Su',
    'Sb', 'Se', 'Si', 'Sn', 'Tb',
    'Te', 'Ur', 'Wt', 'Va', 'Yt', 'Zn',
    'ATM'
]
num2aa.extend(ptms)
num2aa_full.extend(ptms)

to1letter = {
    "ALA":'A', "ARG":'R', "ASN":'N', "ASP":'D', "CYS":'C',
    "GLN":'Q', "GLU":'E', "GLY":'G', "HIS":'H', "ILE":'I',
    "LEU":'L', "LYS":'K', "MET":'M', "PHE":'F', "PRO":'P',
    "SER":'S', "THR":'T', "TRP":'W', "TYR":'Y', "VAL":'V',
    "DA":'a', "DC":'c', "DG":'g', "DT":'t', "DX":'x',
    "RA":'b', "RC":'d', "RG":'h', "RU":'u', "RX":'y',
}
to3letter = {v:k for k,v in to1letter.items()}

aa2num= {x:i for i,x in enumerate(num2aa)}
aa2num_full= {x:i for i,x in enumerate(num2aa_full)}

restype_str_to_int = {
    "A": 0,
    "C": 1,
    "D": 2,
    "E": 3,
    "F": 4,
    "G": 5,
    "H": 6,
    "I": 7,
    "K": 8,
    "L": 9,
    "M": 10,
    "N": 11,
    "P": 12,
    "Q": 13,
    "R": 14,
    "S": 15,
    "T": 16,
    "V": 17,
    "W": 18,
    "Y": 19,
    "X": 20,
}
restype_int_to_str = {
    0: "A",
    1: "C",
    2: "D",
    3: "E",
    4: "F",
    5: "G",
    6: "H",
    7: "I",
    8: "K",
    9: "L",
    10: "M",
    11: "N",
    12: "P",
    13: "Q",
    14: "R",
    15: "S",
    16: "T",
    17: "V",
    18: "W",
    19: "Y",
    20: "X",
}
alphabet = list(restype_str_to_int)
