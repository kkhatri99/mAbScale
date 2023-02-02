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
    Half-body Compositions and Intact Glycoforms.
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built-in Python modules.
from collections import Counter
from itertools import combinations_with_replacement, product
from typing import Dict, List, Tuple

# App modules.
import calc as cal


def calc_halfbody_comp_ctr(hc_pos_1_full_reduced_comp_ctr: Counter,
                           lc_pos_1_full_reduced_comp_ctr: Counter,
                           total_disulfides: int) -> Counter:
    '''
        Calculate the half-body composition.

        Arguments:
            hc_pos_1_full_reduced_comp_ctr {Counter} -- Heavy chain - 1 (clipped / half body) compositon.

            lc_pos_1_full_reduced_comp_ctr {Counter} -- Light chain - 1 (half body) composition.

            total_disulfides {Int} -- Count of total disulfides.

        Returns:
            {Counter} -- Half body composition.
    '''
    molecule_comp_dict: Dict[str, Counter] = cal.gen_molecule_comp_dict()
    hydrogen: Counter = molecule_comp_dict['Hydrogen']

    # Removing all hydrogens from Cys participating in diS bonds.
    halfbody_comp_ctr = hc_pos_1_full_reduced_comp_ctr + \
                        lc_pos_1_full_reduced_comp_ctr - \
                        (cal.mult(int(total_disulfides / 2), cal.mult(2, hydrogen)))

    return halfbody_comp_ctr
    # End of function calc_halfbody_comp_ctr()


def calc_halfbody_glycan_comp_dict(halfbody_comp_ctr: Counter,
                                   halfbody_chain_combo_label: str,
                                   is_lc_glycos: bool,
                                   glycan_comps_dict: Dict[str, Counter]) -> Dict[str, Counter]:
    '''
        Calculate a dictionary of half body chains with glycosylation compositions.

        Arguments:
            halfbody_comp_ctr {Counter} -- Half-body composition.

            halfbody_chain_combo_label {str} -- Label that describes the combination of chains.

            is_lc_glycos {bool} --  If there is light chain glycosolation.

            glycans_comp_dict {Dict} -- Dictionary of glycan compositions.

        Returns:
            {Dict} -- Dictionary of half-body chains with glycosylation compositions.
    '''
    combo_count: int = 0

    combo_count = 2 if is_lc_glycos else 1

    glycan_list: list = list(glycan_comps_dict.keys())
    glycan_combos_list: list = list(combinations_with_replacement(glycan_list, combo_count))

    glycan_combo_name: str = ''
    halfbody_glycan_comp_dict: Dict[str, Counter] = {}
    total_comp_ctr: Counter = Counter({})

    for glycan_combo in glycan_combos_list:
        if combo_count == 2:
            glycan_combo_name = halfbody_chain_combo_label + '|+ ' + \
                                glycan_combo[0] + ' + ' + \
                                glycan_combo[1]
        else:
            glycan_combo_name = halfbody_chain_combo_label + '|+ ' + \
                                glycan_combo[0]

        total_comp_ctr = Counter(halfbody_comp_ctr)

        for cur_index in range(combo_count):
            total_comp_ctr += glycan_comps_dict[glycan_combo[cur_index]]

        halfbody_glycan_comp_dict[glycan_combo_name] = total_comp_ctr

    return halfbody_glycan_comp_dict
    # End of function calc_halfbody_glycan_comp_dict()


def gen_halfbody_glycans_comp_dicts(params_dict: dict,
                                    reduced_comp_dict: dict,
                                    unique_hc_list: list,
                                    unique_lc_list: list,
                                    hc_lc_glycans_comp_dict: dict) -> Tuple[dict, dict, list]:
    '''
        Generate dictionarys for intact reports.

        Arguments:
            params_dict {dict} -- Parameters dictionary.
            reduced_comp_dict {dict} -- Reduced composition properties dictionary.
            unique_hc_list {list} -- Unique heavy chain list.
            unique_lc_list {list} -- Unique light chain list.
            hc_lc_glycans_comp_dict {dict} -- Glycan compositions dictionary for
                                              heavy chain and light chain.

        Returns:
            {Tuple[dict, dict, list]} -- 3 data structures.
                {dict} -- Half-body compositions dictionary.
                {dict} -- Half-body compositions with glycosylation dictionary.
                {list} -- Combinations of heavy chains and light chanins.
    '''
    halfbody_comps_dict: Dict[str, Counter] = {}
    halfbody_glycan_comps_dict: Dict[str, Counter] = {}

    # hc12_lc12_combos: List = cal.hc_lc_combos(hc_chain_list=unique_hc_list,
    #                                           lc_chain_list=unique_lc_list)

    half_body_chains_list: List[List[str]] = [unique_hc_list,
                                              unique_lc_list]
    hc12_lc12_combos: List[Tuple[str, str]] = list(product(*half_body_chains_list))

    for cur_combo in hc12_lc12_combos:
        # halfbody_comps_dict.
        halfbody_combo_label: str = f'{cur_combo[0]}:{cur_combo[1]}'
        halfbody_comp_ctr: Counter = calc_halfbody_comp_ctr(hc_pos_1_full_reduced_comp_ctr=Counter(reduced_comp_dict[cur_combo[0]]['full_reduced_comp_ctr']),
                                                            lc_pos_1_full_reduced_comp_ctr=Counter(reduced_comp_dict[cur_combo[1]]['full_reduced_comp_ctr']),
                                                            total_disulfides=int(params_dict['Total_Disulfides']))
        halfbody_comps_dict[halfbody_combo_label] = halfbody_comp_ctr

        # halfbody_glycan_comps_dict.
        cur_halfbody_glycan_comp_dict: Dict[str, Counter] = calc_halfbody_glycan_comp_dict(halfbody_comp_ctr=halfbody_comp_ctr,
                                                                                           halfbody_chain_combo_label=halfbody_combo_label,
                                                                                           is_lc_glycos=bool(params_dict['LC_Glycos']),
                                                                                           glycan_comps_dict=hc_lc_glycans_comp_dict)
        halfbody_glycan_comps_dict.update(cur_halfbody_glycan_comp_dict)

    return (halfbody_comps_dict, halfbody_glycan_comps_dict, hc12_lc12_combos)
    # End of function gen_halfbody_glycans_comp_dicts
