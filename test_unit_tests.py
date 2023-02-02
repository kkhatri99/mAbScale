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

# pylint: disable=redefined-outer-name
# pylint: disable=missing-function-docstring
'''
    Purpose: Unit tests.
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Build-in modules.
from collections import Counter
from typing import Dict

# Third party modules.
import pytest

# App modules.
import calc as cal
import common_classes as cc
import common_functions as cf
import dialog_box as db


#-------------------------------------------------------------------------------

def test__mult__good():
    new_mol_comp_ctr = cal.mult(multiplier=2,
                              molecule_comp_ctr=Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1}))
    assert new_mol_comp_ctr ==  Counter({'Carbon': 10, 'Hydrogen': 18, 'Nitrogen':2, 'Oxygen': 2, 'Sulfur': 2})


def test__mult__good_with_zeroes():
    new_mol_comp_ctr = cal.mult(multiplier=2,
                              molecule_comp_ctr=Counter({'Carbon': 6, 'Hydrogen': 10, 'Nitrogen': 0, 'Oxygen': 4, 'Sulfur': 0}))
    assert new_mol_comp_ctr ==  Counter({'Carbon': 12, 'Hydrogen': 20, 'Nitrogen':0, 'Oxygen': 8, 'Sulfur': 0})


def test__mult__bad():
    new_mol_comp_ctr = cal.mult(multiplier=2,
                              molecule_comp_ctr=Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1}))
    assert new_mol_comp_ctr !=  Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1})

#---------------------------------------------------------------------------------------------------

@pytest.fixture
def tst_ele_mass_dict():
    ele_name_mass_dict: Dict[str, Dict[str, float]] = {}
    ele_name_mass_dict = {'Carbon': {'Mass': 1.0},
                          'Hydrogen': {'Mass': 2.0},
                          'Nitrogen': {'Mass': 3.0},
                          'Oxygen': {'Mass': 4.0},
                          'Sulfur': {'Mass': 5.1}}

    return ele_name_mass_dict


def test__calc_mass_by_ele_name__good(tst_ele_mass_dict):
    new_mass = cf.calc_mass_by_ele_name(mol_comp_ctr=Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1}),
                                        ele_name_mass_dict=tst_ele_mass_dict)
    assert new_mass == 35.1


def test__calc_mass_by_ele_name__good_with_zeroes(tst_ele_mass_dict):
    new_mass = cf.calc_mass_by_ele_name(mol_comp_ctr=Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':0, 'Oxygen': 1, 'Sulfur': 1}),
                                        ele_name_mass_dict=tst_ele_mass_dict)
    assert new_mass == 32.1


def test__calc_mass_by_ele_name__bad(tst_ele_mass_dict):
    new_mass = cf.calc_mass_by_ele_name(mol_comp_ctr=Counter({'Carbon': 5, 'Hydrogen': 9, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1}),
                                        ele_name_mass_dict=tst_ele_mass_dict)
    assert new_mass != 17.2

# ------------------------------------------------------------------------------

def test__are_chars_in_list__good():
    result = db.are_chars_in_list_valid(check_string='ABC',
                                        valid_list=['A', 'B', 'C', 'X', 'Y', 'Z'])
    assert result is True


def test__are_chars_in_list__good_str_case():
    result = db.are_chars_in_list_valid(check_string='AbC',
                                        valid_list=['A', 'B', 'C', 'X', 'Y', 'Z'])
    assert result is True


def test__are_chars_in_list__good_list_case():
    result = db.are_chars_in_list_valid(check_string='ABC',
                                        valid_list=['A', 'B', 'c', 'X', 'Y', 'Z'])
    assert result is True


def test__are_chars_in_list__bad():
    result = db.are_chars_in_list_valid(check_string='ABD',
                                        valid_list=['A', 'B', 'C', 'X', 'Y', 'Z'])
    assert result is False

# ------------------------------------------------------------------------------

def test__get_seq_elements__good():
    new_mol_comp_ctr = cal.get_seq_ele_comp(sequence='PEP')
    assert new_mol_comp_ctr == Counter({'Carbon': 15, 'Hydrogen': 21, 'Nitrogen':3, 'Oxygen': 5})


def test__get_seq_elements__good_case():
    new_mol_comp_ctr = cal.get_seq_ele_comp(sequence='tides')
    assert new_mol_comp_ctr == Counter({'Carbon': 22, 'Hydrogen': 35, 'Nitrogen':5, 'Oxygen': 11})


def test__get_seq_elements__bad():
    new_mol_comp_ctr = cal.get_seq_ele_comp(sequence='PEP')
    assert new_mol_comp_ctr != Counter({'Carbon': 13, 'Hydrogen': 29, 'Nitrogen':3, 'Oxygen': 5})

# ------------------------------------------------------------------------------

def test__gen_amino_acid_dict__good_count():
    aa_dict = cf.gen_amino_acid_dict()
    aa_count = len(aa_dict)
    assert aa_count == 20


def test__gen_amino_acid_dict__good_specific_key():
    aa_dict = cf.gen_amino_acid_dict()
    current_aa = aa_dict.get('C')
    assert current_aa == 'Cys'


def test__gen_amino_acid_dict__bad_invalid_key():
    aa_dict = cf.gen_amino_acid_dict()
    current_aa = aa_dict.get('Rhino', '_invalid_')
    assert current_aa == '_invalid_'


def test__gen_amino_acid_dict__bad_specific_key():
    aa_dict = cf.gen_amino_acid_dict()
    current_aa = aa_dict.get('A')
    assert current_aa != 'Dog'

# ------------------------------------------------------------------------------

def test__load_ele_mass_dict__good_count():
    ele_name_mass_dict = cf.gen_ele_name_mass_dict()
    ele_mass_count = len(ele_name_mass_dict)
    assert ele_mass_count == 15


def test__load_ele_mass_dict__good_specific_key_1():
    ele_name_mass_dict = cf.gen_ele_name_mass_dict()
    current_ele = ele_name_mass_dict.get('Carbon', {'Mass': -999})
    current_ele_mass = current_ele.get('Mass')
    assert current_ele_mass == 12.0107359


def test__load_ele_mass_dict__good_specific_key_2():
    ele_name_mass_dict = cf.gen_ele_name_mass_dict()
    current_ele_mass = ele_name_mass_dict['Oxygen']['Mass']
    assert current_ele_mass == 15.99940492


def test__load_ele_mass_dict__bad_invalid_key():
    ele_name_mass_dict = cf.gen_ele_name_mass_dict()
    current_ele_mass = ele_name_mass_dict.get('Rhino', '_invalid_')
    assert current_ele_mass == '_invalid_'


def test__load_ele_mass_dict__bad_specific_key():
    ele_name_mass_dict = cf.gen_ele_name_mass_dict()
    current_ele_mass = ele_name_mass_dict['Hydrogen']['Mass']
    assert current_ele_mass != -1.2345

# ------------------------------------------------------------------------------

def test__load_glycan_dict__good_count():
    glycan_comp_dict = cal.load_glycan_dict(chain_type=cc.Chain.HC)
    glycan_count = len(glycan_comp_dict)
    assert glycan_count == 6


def test__load_glycan_dict__good_specific_key():
    glycan_comp_dict = cal.load_glycan_dict(chain_type=cc.Chain.HC)
    current_glycan = glycan_comp_dict.get('G0F')
    assert current_glycan == Counter({'Hex': 3, 'HexNAc': 4, 'NeuAc': 0, 'NeuGc': 0, 'dHex': 1})


def test__load_glycan_dict__bad_invalid_key():
    glycan_comp_dict = cal.load_glycan_dict(chain_type=cc.Chain.HC)
    current_glycan = glycan_comp_dict.get('Rhino', '_invalid_')
    assert current_glycan == '_invalid_'


def test__load_glycan_dict__bad_specific_key():
    glycan_comp_dict = cal.load_glycan_dict(chain_type=cc.Chain.HC)
    current_glycan = glycan_comp_dict.get('G0')
    assert current_glycan != Counter({'Hex': 55, 'HexNAc': 44, 'NeuAc': 33, 'NeuGc': 22, 'dHex': 11})

# ------------------------------------------------------------------------------

def test__gen_glycan_comp_dict__good_count():
    glycan_comp_dict = cal.gen_glycan_comp_dict(chain_type=cc.Chain.HC)
    glycan_count = len(glycan_comp_dict)
    assert glycan_count == 6


def test__gen_glycan_comp_dict__good_specific_key():
    glycan_comp_dict = cal.gen_glycan_comp_dict(chain_type=cc.Chain.HC)
    current_glycan = glycan_comp_dict.get('G0F')
    assert current_glycan == Counter({'Carbon': 56, 'Hydrogen': 92, 'Nitrogen': 4, 'Oxygen': 39})


def test__gen_glycan_comp_dict__bad_invalid_key():
    glycan_comp_dict = cal.gen_glycan_comp_dict(chain_type=cc.Chain.HC)
    current_glycan = glycan_comp_dict.get('Rhino', '_invalid_')
    assert current_glycan == '_invalid_'


def test__gen_glycan_comp_dict__bad_specific_key():
    glycan_comp_dict = cal.gen_glycan_comp_dict(chain_type=cc.Chain.HC)
    current_glycan = glycan_comp_dict.get('G0')
    assert current_glycan != Counter({'Hydrogen': 1, 'Carbon': 1, 'Oxygen': 1, 'Nitrogen': 1})

# ------------------------------------------------------------------------------

def test__gen_molecule_comp_dict__count():
    mol_comp_dict = cal.gen_molecule_comp_dict()
    mol_ele_count = len(mol_comp_dict)
    assert mol_ele_count == 31


def test__gen_molecule_comp_dict__good_specific_key():
    mol_comp_dict = cal.gen_molecule_comp_dict()
    current_mol_comp_ctr = mol_comp_dict.get('Gln')
    assert current_mol_comp_ctr == Counter({'Carbon': 5, 'Hydrogen': 8, 'Nitrogen':2, 'Oxygen': 2, 'Sulfur': 0})


def test__gen_molecule_comp_dict__bad_invalid_key():
    mol_comp_dict = cal.gen_molecule_comp_dict()
    current_mol_comp_ctr = mol_comp_dict.get('Rhino', '_invalid_')
    assert current_mol_comp_ctr == '_invalid_'


def test__gen_molecule_comp_dict__bad_specific_key():
    mol_comp_dict = cal.gen_molecule_comp_dict()
    current_mol_comp_ctr = mol_comp_dict.get('Gln')
    assert current_mol_comp_ctr != Counter({'Carbon': 1, 'Hydrogen': 1, 'Nitrogen':1, 'Oxygen': 1, 'Sulfur': 1})

# ------------------------------------------------------------------------------

def test__is_int__good_str():
    assert db.is_int(value='2') is True


def test__is_int__bad_str():
    assert db.is_int(value='cat') is False


def test__is_int__bad_type():
    assert db.is_int(value=str(3.3)) is False

# ------------------------------------------------------------------------------

def test__is_even__good_str():
    assert db.is_even(num=2) is True


def test__is_even__bad_str():
    assert db.is_even(num=7) is False


def test__is_even__good_zero():
    assert db.is_even(num=0) is True

# ------------------------------------------------------------------------------

def test__is_in_range__good():
    assert db.is_in_range_int(num=5,
                              min_val=1,
                              max_val=10) is True

def test__is_in_range__good_edge_min():
    assert db.is_in_range_int(num=1,
                              min_val=1,
                              max_val=10) is True

def test__is_in_range__good_edge_max():
    assert db.is_in_range_int(num=10,
                              min_val=1,
                              max_val=10) is True

def test__is_in_range__bad_edge_min():
    assert db.is_in_range_int(num=0,
                              min_val=1,
                              max_val=10) is False

def test__is_in_range__bad_edge_max():
    assert db.is_in_range_int(num=11,
                              min_val=1,
                              max_val=10) is False

# ------------------------------------------------------------------------------
