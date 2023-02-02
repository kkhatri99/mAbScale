"""
Copyright 2022 GlaxoSmithKline Research & Development Limited

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
"""

'''
    Common functions
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built-in modules.
from collections import Counter
import os
import re
import sys
from typing import Dict

# Third-party modules.
import pandas as pd


def gen_amino_acid_dict() -> Dict[str, str]:
    '''
        Generate a dictionary of Amino Acids.

        Used for conversion of 1 letter abbreviation to 3 letter abbreviation.

        Arguments:
            None

        Returns:
            {Dict} -- Dictionary of Amino Acids.

        Reference:
            https://www.ddbj.nig.ac.jp/ddbj/code-e.html
    '''
    aa1_aa3_dict: Dict[str, str] = {}

    aa1_aa3_dict['A'] = 'Ala'
    aa1_aa3_dict['R'] = 'Arg'
    aa1_aa3_dict['N'] = 'Asn'
    aa1_aa3_dict['D'] = 'Asp'
    aa1_aa3_dict['C'] = 'Cys'
    aa1_aa3_dict['Q'] = 'Gln'
    aa1_aa3_dict['E'] = 'Glu'
    aa1_aa3_dict['G'] = 'Gly'
    aa1_aa3_dict['H'] = 'His'
    aa1_aa3_dict['I'] = 'Ile'
    aa1_aa3_dict['L'] = 'Leu'
    aa1_aa3_dict['K'] = 'Lys'
    aa1_aa3_dict['M'] = 'Met'
    aa1_aa3_dict['F'] = 'Phe'
    aa1_aa3_dict['P'] = 'Pro'
    # aa1_aa3_dict['O'] = 'Pyl'
    aa1_aa3_dict['S'] = 'Ser'
    # aa1_aa3_dict['U'] = 'Sec'
    aa1_aa3_dict['T'] = 'Thr'
    aa1_aa3_dict['W'] = 'Trp'
    aa1_aa3_dict['Y'] = 'Tyr'
    aa1_aa3_dict['V'] = 'Val'
    # aa1_aa3_dict['B'] = 'Asx'
    # aa1_aa3_dict['Z'] = 'Glx'
    # aa1_aa3_dict['X'] = 'Xaa'
    # aa1_aa3_dict['J'] = 'Xle'

    return aa1_aa3_dict
    # End of function gen_amino_acid_dict()


def formula_to_composition(mol_formula: str) -> Counter:
    '''
        Converts a molecular formula to it's elemental composition.

        Arguments:
            mol_formula {str} -- Molecular formula
            ex: 'H2O'

        Returns:
            {Counter} -- Elemental composition.
            ex: Counter({'H': 2, 'O': 1})

        Notes:
            1. element_name is: capital letter followed by optional lower-case letter.
            2. count is: empty string (so the count is 1), or a set of digits
    '''
    element_pat = re.compile(R"([A-Z][a-z]?)(\d*)")
    all_elements = []

    for (element_name, count) in element_pat.findall(mol_formula):
        count = 1 if count == '' else int(count)
        all_elements.extend([element_name] * count)

    molecule_comp_ctr = Counter(all_elements)
    return molecule_comp_ctr
    # End of function formula_to_composition()


def load_ele_mass_df() -> pd.DataFrame:
    '''
        Loads element masses information from the "Element_Mass.csv" file into
        a pandas DataFrame.

        Arguments:
            None

        Returns:
            {DataFrame} -- DataFrame of element masses.
    '''
    application_path: str = ''

    if getattr(sys, 'frozen', False):
        application_path = os.path.dirname(sys.executable)
    elif __file__:
        application_path = os.path.dirname(__file__)

    file_name_path: str = application_path + R'\Element_Mass.csv'
    ele_mass_df = pd.read_csv(file_name_path)
    return ele_mass_df
    # End of function load_ele_mass_df()


def gen_ele_symbol_mass_dict() -> Dict[str, Dict[str, float]]:
    '''
        Converts element mass DataFrame into dictionary with element symbol as key.

        Arguments:
            None

        Returns:
            {Dict} -- Dictionary of element masses.

        Example Use:
            In:

            Out:
                {
                    'Br': {'Mass': 79.904},
                    'Ca': {'Mass': 39.9625912},
                    'C': {'Mass': 12.0107359},
                    ...
                }
    '''
    ele_mass_df: pd.DataFrame = load_ele_mass_df()
    ele_symbol_mass_dict = ele_mass_df.set_index('Symbol')[['Mass']].T.to_dict('dict')
    return ele_symbol_mass_dict
    # End of function gen_ele_symbol_mass_dict()


def gen_ele_name_mass_dict() -> Dict[str, Dict[str, float]]:
    '''
        Converts element mass DataFrame into dictionary with element name as key.

        Arguments:
            None

        Returns:
            {Dict} -- Dictionary of element masses.

        Example Use:
            In:

            Out:
                {
                    'Bromine': {'Mass': 79.904},
                    'Calcium': {'Mass': 39.9625912},
                    'Carbon': {'Mass': 12.0107359},
                    ...
                }
    '''
    ele_mass_df: pd.DataFrame = load_ele_mass_df()
    ele_name_mass_dict = ele_mass_df.set_index('Name')[['Mass']].T.to_dict('dict')
    return ele_name_mass_dict
    # End of function gen_ele_name_mass_dict()


def calc_mass_by_ele_symbol(mol_comp_ctr: Counter,
                            ele_symbol_mass_dict: Dict[str, Dict[str, float]]) -> float:
    '''
        Calculate mass of a molecule based by element symbol(s).

        Arguments:
            mol_comp_ctr {Counter} -- Molecule to calculate mass of.

            ele_symbol_mass_dict {Dict[str, float]} -- Dictionary of element masses.

        Returns:
            {float} -- Mass of molecule.
    '''
    mass: float = 0.0

    for key, value in mol_comp_ctr.items():
        mass += ele_symbol_mass_dict[key]['Mass'] * value

    return mass
    # End of function calc_mass_by_ele_symbol()


def calc_mass_by_ele_name(mol_comp_ctr: Counter,
                          ele_name_mass_dict: Dict[str, Dict[str, float]]) -> float:
    '''
        Calculate mass of a given molecule by element name(s).

        Arguments:
            mol_comp_ctr {Counter} -- Molecule to calculate mass of.

            ele_name_mass_dict {Dict[str, float]} -- Dictionary of element masses.

        Returns:
            {float} -- Mass of molecule.
    '''
    mass: float = 0.0

    for key, value in mol_comp_ctr.items():
        mass += ele_name_mass_dict[key]['Mass'] * value

    return mass
    # End of function calc_mass_by_ele_name()
