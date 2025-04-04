# Copyright 2021 AlQuraishi Laboratory
# Copyright 2021 DeepMind Technologies Limited
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#      http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

# This mapping is used when we need to store atom data in a format that requires
# fixed atom data size for every residue (e.g. a numpy array).
atom_types = [
    "N",
    "CA",
    "C",
    "CB",
    "O",
    "CG",
    "CG1",
    "CG2",
    "OG",
    "OG1",
    "SG",
    "CD",
    "CD1",
    "CD2",
    "ND1",
    "ND2",
    "OD1",
    "OD2",
    "SD",
    "CE",
    "CE1",
    "CE2",
    "CE3",
    "NE",
    "NE1",
    "NE2",
    "OE1",
    "OE2",
    "CH2",
    "NH1",
    "NH2",
    "OH",
    "CZ",
    "CZ2",
    "CZ3",
    "NZ",
    "OXT",
]
atom_order = {atom_type: i for i, atom_type in enumerate(atom_types)}
atom_type_num = len(atom_types)  # := 37.

restype_1to3 = {
    "A": "ALA",
    "R": "ARG",
    "N": "ASN",
    "D": "ASP",
    "C": "CYS",
    "Q": "GLN",
    "E": "GLU",
    "G": "GLY",
    "H": "HIS",
    "I": "ILE",
    "L": "LEU",
    "K": "LYS",
    "M": "MET",
    "F": "PHE",
    "P": "PRO",
    "S": "SER",
    "T": "THR",
    "W": "TRP",
    "Y": "TYR",
    "V": "VAL",
}

# Approximate Volumes of amino acids in cubic angstroms.
# https://www.imgt.org/IMGTeducation/Aide-memoire/_UK/aminoacids/abbreviation.html
amino_acid_volumes = {
    "A": 88.6,  # Alanine
    "R": 173.4,  # Arginine
    "N": 114.1,  # Asparagine
    "D": 111.1,  # Aspartic acid
    "C": 108.5,  # Cysteine
    "Q": 143.8,  # Glutamine
    "E": 138.4,  # Glutamic acid
    "G": 60.1,  # Glycine
    "H": 153.2,  # Histidine
    "I": 166.7,  # Isoleucine
    "L": 166.7,  # Leucine
    "K": 168.6,  # Lysine
    "M": 162.9,  # Methionine
    "F": 189.9,  # Phenylalanine
    "P": 112.7,  # Proline
    "S": 89.0,  # Serine
    "T": 116.1,  # Threonine
    "W": 227.8,  # Tryptophan
    "Y": 193.6,  # Tyrosine
    "V": 140.0,  # Valine
    "X": 88.6,  # Unknown, use Alanine as approximation
}
