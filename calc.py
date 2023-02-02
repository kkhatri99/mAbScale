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
    Functions and classes for performing calculations.
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built-in modules.
from collections import Counter
from itertools import combinations_with_replacement
import os
import sys
from typing import Dict, List

# Third party modules.
import pandas as pd

# App modules.
import common_classes as cc
import common_functions as cf


def gen_molecule_comp_dict() -> Dict[str, Counter]:
    '''
        Generate a dictionary of molecule compositions.

        Arguments:
            None

        Returns:
            {Dict} -- Dictionary of molecule compositions.
    '''
    mol_comp_dict: Dict[str, Counter] = {}

    # Elements
    mol_comp_dict['Carbon'] = Counter({'Carbon': 1, 'Hydrogen': 0, 'Nitrogen': 0, 'Oxygen': 0, 'Sulfur': 0})
    mol_comp_dict['Hydrogen'] = Counter({'Carbon': 0, 'Hydrogen': 1, 'Nitrogen': 0, 'Oxygen': 0, 'Sulfur': 0})
    mol_comp_dict['Nitrogen'] = Counter({'Carbon': 0, 'Hydrogen': 0, 'Nitrogen': 1, 'Oxygen': 0, 'Sulfur': 0})
    mol_comp_dict['Oxygen'] = Counter({'Carbon': 0, 'Hydrogen': 0, 'Nitrogen': 0, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Sulfur'] = Counter({'Carbon': 0, 'Hydrogen': 0, 'Nitrogen': 0, 'Oxygen': 0, 'Sulfur': 1})

    # Misc.
    mol_comp_dict['Water'] = Counter({'Carbon':0,'Hydrogen':2,'Nitrogen':0,'Oxygen':1,'Sulfur':0})

    # Amino Acid Residue - Residue Formula.
    mol_comp_dict['Ala'] = Counter({'Carbon': 3, 'Hydrogen': 5, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Arg'] = Counter({'Carbon': 6, 'Hydrogen': 12, 'Nitrogen':4, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Asn'] = Counter({'Carbon': 4, 'Hydrogen': 6, 'Nitrogen':2, 'Oxygen': 2, 'Sulfur': 0})
    mol_comp_dict['Asp'] = Counter({'Carbon': 4, 'Hydrogen': 5, 'Nitrogen':1, 'Oxygen': 3, 'Sulfur': 0})
    mol_comp_dict['Cys'] = Counter({'Carbon': 3, 'Hydrogen': 5, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1})
    mol_comp_dict['Glu'] = Counter({'Carbon': 5, 'Hydrogen': 7, 'Nitrogen':1, 'Oxygen': 3, 'Sulfur': 0})
    mol_comp_dict['Gln'] = Counter({'Carbon': 5, 'Hydrogen': 8, 'Nitrogen':2, 'Oxygen': 2, 'Sulfur': 0})
    mol_comp_dict['Gly'] = Counter({'Carbon': 2, 'Hydrogen': 3, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['His'] = Counter({'Carbon': 6, 'Hydrogen': 7, 'Nitrogen':3, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Ile'] = Counter({'Carbon': 6, 'Hydrogen': 11, 'Nitrogen':1,'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Leu'] = Counter({'Carbon': 6, 'Hydrogen': 11, 'Nitrogen':1,'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Lys'] = Counter({'Carbon': 6, 'Hydrogen': 12, 'Nitrogen':2, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Met'] = Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1})
    mol_comp_dict['Phe'] = Counter({'Carbon': 9, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Pro'] = Counter({'Carbon': 5, 'Hydrogen': 7, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 0})
    mol_comp_dict['Ser'] = Counter({'Carbon': 3, 'Hydrogen': 5, 'Nitrogen':1, 'Oxygen': 2, 'Sulfur': 0})
    mol_comp_dict['Thr'] = Counter({'Carbon': 4, 'Hydrogen': 7, 'Nitrogen':1, 'Oxygen': 2, 'Sulfur': 0})
    mol_comp_dict['Trp'] = Counter({'Carbon': 11,'Hydrogen': 10, 'Nitrogen':2, 'Oxygen': 1,'Sulfur': 0})
    mol_comp_dict['Tyr'] = Counter({'Carbon': 9, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 2, 'Sulfur': 0})
    mol_comp_dict['Val'] = Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 0})

    # Monosacarides - Residue Formula.
    mol_comp_dict['HexNAc'] =  Counter({'Carbon': 8, 'Hydrogen': 13, 'Nitrogen': 1, 'Oxygen': 5, 'Sulfur': 0})
    mol_comp_dict['Hex'] = Counter({'Carbon': 6, 'Hydrogen': 10, 'Nitrogen': 0, 'Oxygen': 5, 'Sulfur': 0})
    mol_comp_dict['dHex'] =  Counter({'Carbon': 6, 'Hydrogen': 10, 'Nitrogen': 0, 'Oxygen': 4, 'Sulfur': 0})
    mol_comp_dict['NeuAc'] = Counter({'Carbon': 11, 'Hydrogen': 17, 'Nitrogen': 1, 'Oxygen': 8, 'Sulfur': 0})
    mol_comp_dict['NeuGc'] = Counter({'Carbon': 11, 'Hydrogen': 17, 'Nitrogen': 1, 'Oxygen': 9, 'Sulfur': 0})

    return mol_comp_dict
    # End of function gen_molecule_comp_dict()


def get_seq_ele_comp(sequence: str) -> Counter:
    '''
        Gets the element composition for the given sequence.

        Arguments:
            sequence {str} -- Peptide sequence.

        Returns:
            {Counter} -- Molecule composition.
    '''
    aa1_aa3_dict: Dict[str, str] = cf.gen_amino_acid_dict()
    molecule_comp_dict: Dict[str, Counter] = gen_molecule_comp_dict()

    sequence = sequence.strip().upper()

    seq_comp_ctr: Counter = Counter({})

    for aa_abbrev1 in sequence:
        aa_abbrev3: str = aa1_aa3_dict[aa_abbrev1]
        aa_comp_ctr: Counter = molecule_comp_dict[aa_abbrev3]
        seq_comp_ctr += aa_comp_ctr

    return seq_comp_ctr
    # End of fucntion get_seq_ele_comp()


def mult(multiplier: int,
         molecule_comp_ctr: Counter) -> Counter:
    '''
        Multiple the number of elements in a molecule composition by a given number.

        Arguments:
            multiplier {int} -- Number to mutiply by.

            exist_dict {Counter} -- A molecule composition.

        Returns:
            {Counter} -- New molecule composition.
    '''
    return Counter({key: value * multiplier for key, value in molecule_comp_ctr.items()})
    # End of function mult()


def load_glycan_dict(chain_type: str) -> Dict:
    '''
        Loads glycan information from "Glycans.csv" file into a dictionary.
        Glycans are filtered and sorted.

        Arguments:
            chain_type {cc.Chain} -- Chain type of HC or LC.
                ex: cc.Chain.HC

        Returns:
            {Dict} -- Dictionary of glycan compositions.
    '''
    application_path: str = ''

    if getattr(sys, 'frozen', False):
        application_path = os.path.dirname(sys.executable)
    elif __file__:
        application_path = os.path.dirname(__file__)

    file_name_path: str = application_path + R'\Glycans.csv'
    glycan_dict: dict = {}

    glycans_df = pd.read_csv(file_name_path)
    if chain_type == cc.Chain.HC:
        filt_sort_df = glycans_df[(glycans_df['Show_HC'].str.upper() == 'Y')].sort_values(by=['Ord_Pos'], ascending=True)
    else:
        filt_sort_df = glycans_df[(glycans_df['Show_LC'].str.upper() == 'Y')].sort_values(by=['Ord_Pos'], ascending=True)

    glycan_dict = filt_sort_df.set_index('Glycan_Name')[['HexNAc', 'Hex', 'dHex', 'NeuAc', 'NeuGc']].T.to_dict('dict')

    return glycan_dict
    # End of function load_glycan_dict()


def gen_glycan_comp_dict(chain_type: str) -> Dict[str, Counter]:
    '''
        Generates glycan composition dictionary.

        Arguments:
            chain_type {str} -- Chain type of HC or LC.

        Returns:
            {Dict} -- Dictionary of glycan compositions.
    '''
    molecule_comp_dict: Dict[str, Counter] = gen_molecule_comp_dict()
    hex_nac = molecule_comp_dict['HexNAc']
    hex_ = molecule_comp_dict['Hex']
    d_hex = molecule_comp_dict['dHex']
    neu_ac = molecule_comp_dict['NeuAc']
    neu_gc = molecule_comp_dict['NeuGc']

    glycan_dict = load_glycan_dict(chain_type=chain_type)
    glycan_comp_dict: Dict[str, Counter] = {}

    for key, value in glycan_dict.items():
        glycan_comp_dict[key] = mult(value['HexNAc'], hex_nac) + \
                                mult(value['Hex'], hex_) + \
                                mult(value['dHex'], d_hex) + \
                                mult(value['NeuAc'], neu_ac) + \
                                mult(value['NeuGc'], neu_gc)

    return glycan_comp_dict
    # End of fucntion gen_glycan_comp_dict()


class HCLCCombo:
    '''
        Container for heavy chain and light chain combo.
    '''
    def __init__(self):
        self.hc_pos_1: int
        self.hc_pos_2: int
        self.lc_pos_1: int
        self.lc_pos_2: int
    # End of class HCLCCombo


def hc_lc_combos(hc_chain_list: list,
                 lc_chain_list: list) -> list:
    '''
        Determines heavy chain and light chain combinations.

        Arguments:
            hc_chain_list {list} -- List of heavy chains.
                ex: ['HC-1', 'HC-2']

            lc_chain_list {list} -- List of light chains.
                ex: ['LC-1']

        Returns:
            {list} -- Heavy chain and light chain combinations.
                ex:  [['HC-1', 'HC-1', 'LC-1', 'LC-1'],
                      ['HC-1', 'HC-2', 'LC-1', 'LC-1'],
                      ...
                     ]
    '''
    hc_combo_list: list = []
    lc_combo_list: list = []

    hc_chain_combos_list: list = list(combinations_with_replacement(hc_chain_list, 2))
    for cur_chain_combo in hc_chain_combos_list:
        hc_combo_list.append(','.join(cur_chain_combo))

    lc_chain_combos_list: list = list(combinations_with_replacement(lc_chain_list, 2))
    for cur_chain_combo in lc_chain_combos_list:
        lc_combo_list.append(','.join(cur_chain_combo))

    hc_lc_combos_list = [f'{hc},{lc}' for hc in hc_combo_list for lc in lc_combo_list]

    combos: list = []

    for cur_combo in hc_lc_combos_list:
        chain_combo = cur_combo.split(',')
        combos.append(chain_combo)

    return combos
    # End of function hc_lc_combos()


def build_chem_mod_dict(chains_list: List,
                        params_dict: Dict) -> Dict:
    '''
        Builds a chem modificaton dictionary wich contains formula, composition,
        and mass.

        Parameters:
            chains_list {list} -- List of chains.
                ex: ['HC-1', 'HC-2', 'LC-1', 'LC-2']

            params_dict {dict} -- User input parameters.

        Returns:
            {Dict} - Dictionary containing formual, composition, and mass for
                     each chain.
            ex:     'HC-1': {'Formula': 'H2O',
                             'Comp': Counter({'H': 2, 'O': 1}),
                             'Mass': 18.01528642},
                    ...
    '''
    ele_symbol_mass_dict = cf.gen_ele_symbol_mass_dict()
    chem_mod_dict: dict = {}

    for chain in chains_list:
        chem_mod_formula: str = params_dict[f'{chain}_Chem_Mod']
        chem_mod_comp_ctr: Counter = cf.formula_to_composition(chem_mod_formula)
        chem_mod_mass: float = cf.calc_mass_by_ele_symbol(mol_comp_ctr=chem_mod_comp_ctr,
                                                          ele_symbol_mass_dict=ele_symbol_mass_dict)
        
        if len(chem_mod_formula) > 0 and chem_mod_formula[0] == "-":
            chem_mod_mass *= -1

        chem_mod_dict[chain] = {'Formula': chem_mod_formula,
                                'Comp': chem_mod_comp_ctr,
                                'Mass': chem_mod_mass}

    return chem_mod_dict
    # End of function build_chem_mod_dict()


def build_chem_mod_chain_combo_dict(chain_combo_list: list,
                                    num_of_chains: int,
                                    chem_mod_dict: dict) -> Dict:
    '''
        Creates chemical modification dictionary for chain combinations - Formual and Mass.

        Aurguemnts:
            combo_list {list} -- list of list of chain combinations.
                Ex: [['HC-1', 'HC-1', 'LC-1', 'LC-2'],
                     ['HC-1', 'HC-2', 'LC-1', 'LC-2'],...]

            num_of_chains {int} -- Number of chains in a combinations.
                Note: The only valid values are 2 and 4.
                Ex: 4

            chem_mod_dict {dict} -- Dictionary of chemical modifications' formula and mass.

        Returns:
            {Dict}
    '''
    chem_mod_chains_combo_dict: dict = {}

    # Check the num_of_chains parameter.
    if num_of_chains not in [2, 4]:
        return chem_mod_chains_combo_dict

    for chain_combo in chain_combo_list:
        chain_combo_len: int = len(chain_combo)

        if chain_combo_len == 4:
            label: str = f'{chain_combo[0]}/{chain_combo[1]}:{chain_combo[2]}/{chain_combo[3]}'
        elif chain_combo_len == 2:
            label: str = f'{chain_combo[0]}:{chain_combo[1]}'
        else: # Bad combination of chains.
            continue

        chem_mod_text: str = ''
        total_chem_mod_mass: float = 0.0

        for chain in chain_combo:
            chem_mod_text += chem_mod_dict[chain]['Formula'] + '; '
            total_chem_mod_mass += chem_mod_dict[chain]['Mass']

        chem_mod_chains_combo_dict[label] = {'Formulas': chem_mod_text[:-2],
                                             'Mass': total_chem_mod_mass}

    return chem_mod_chains_combo_dict
    # End of function build_chem_mod_chain_combo_dict()


def build_chain_dict(chains_list: list,
                     params_dict: dict) -> Dict:
    '''
        Builds a chain dictionary with suequence and compostion properties.

        Parameters:
            chains_list {list} -- List of chains.
                ex: ['HC-1', 'HC-2', 'LC-1', 'LC-2']

            params_dict {dict} -- User input parameters.

        Returns:
            {Dict} -- Dictionary of each chains sequence and composition.
                ex:  'HC-1': {'Sequence': 'PEP',
                              'Comp': Counter({'Hydrogen': 23, 'Carbon': 15, 'Oxygen': 6, 'Nitrogen': 3})}, ...
    '''
    chain_dict: dict = {}

    molecule_comp_dict: Dict[str, Counter] = gen_molecule_comp_dict()
    water: Counter = molecule_comp_dict['Water']

    for chain in chains_list:
        sequence: str = params_dict[f'{chain}_Sequence']
        comp_ctr: Counter = get_seq_ele_comp(sequence=sequence)
        comp_ctr += water
        chain_dict[chain] = {'Sequence': sequence,
                              'Comp': comp_ctr}

    return chain_dict
    # End of function build_chain_dict()


def gen_unique_chain_list(chain_type: str,
                          params_dict: Dict) -> List[str]:
    '''
        Produces a list of unique chains based on sequence.

        Aurguments:
            chain_type {str} -- Type of chain - HC or LC.
                ex: cc.Chain.HC

            params_dict {Dict} -- Parameters dictionary.


        Returns
            {List[str]} -- List of unqiue chains.
                ex 1: ['HC-1', 'HC-2']
                ex 2: ['LC-1']
    '''
    unique_chain_list: List[str] = []

    unique_chain_list.append(f'{chain_type}-1')

    if params_dict[f'{chain_type}-2_Sequence'] != params_dict[f'{chain_type}-1_Sequence']:
        unique_chain_list.append(f'{chain_type}-2')

    return unique_chain_list
    # End of function gen_unique_chain_list()
