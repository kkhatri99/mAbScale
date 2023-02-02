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
    Intact Compositions and Intact Glycoforms.
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built-in modules.
from collections import Counter
from itertools import combinations_with_replacement
from typing import Dict, List, Tuple

# App modules.
import calc as cal


def calc_intact_comp_ctr(hc_pos_1_full_reduced_comp_ctr: Counter,
                         hc_pos_2_full_reduced_comp_ctr: Counter,
                         lc_pos_1_full_reduced_comp_ctr: Counter,
                         lc_pos_2_full_reduced_comp_ctr: Counter,
                         total_disulfides: int) -> Counter:
    '''
        Calculate the intact composition.

        Arguments:
            hc_pos_1_full_reduced_comp_ctr {Counter} -- Heavy chain - 1 (clipped / intact) compositon.

            hc_pos_2_full_reduced_comp_ctr {Counter} -- Heavy chain -2 (clipped / intact) composition.

            lc_pos_1_full_reduced_comp_ctr {Counter} -- Light chain - 1 (intact) composition.

            lc_pos_2_full_reduced_comp_ctr {Counter} -- Light chain - 2 (intact) composition.

            total_disulfides {Int} -- Count of total disulfides.

        Returns:
            {Counter} -- Intact composition.
    '''
    molecule_comp_dict: Dict[str, Counter] = cal.gen_molecule_comp_dict()
    hydrogen: Counter = molecule_comp_dict['Hydrogen']

    # Removing all hydrogens from Cys participating in diS bonds.
    intact_comp_ctr = hc_pos_1_full_reduced_comp_ctr + \
                      hc_pos_2_full_reduced_comp_ctr + \
                      lc_pos_1_full_reduced_comp_ctr + \
                      lc_pos_2_full_reduced_comp_ctr - \
                      (cal.mult(int(total_disulfides), cal.mult(2, hydrogen)))

    return intact_comp_ctr
    # End of function calc_intact_comp_ctr()


def calc_intact_glycan_comp_dict(intact_comp_ctr: Counter,
                                 intact_chain_combo_label: str,
                                 is_lc_glycos: bool,
                                 glycan_comps_dict: Dict[str, Counter]) -> Dict[str, Counter]:
    '''
        Calculate a dictionary of intact chain with glycan compositions.

        Arguments:
            intact_comp_ctr {Counter} -- Intact composition.

            intact_chain_combo_label {str} -- Label that describes the combination of chains.

            is_lc_glycos {bool} -- If there is light chain glycosolation.

            glycans_comp_dict {Dict} -- Dictionary of glycan compositions.

        Returns:
            {Dict} -- Dictionary of intact chain with glycans compositions.
    '''
    combo_count: int = 0

    combo_count = 4 if is_lc_glycos else 2

    glycan_list: list = list(glycan_comps_dict.keys())
    glycan_combos_list: list = list(combinations_with_replacement(glycan_list, combo_count))

    glycan_combo_name: str = ''
    intact_glycan_comp_dict: Dict[str, Counter] = {}
    total_comp_ctr: Counter = Counter({})

    for glycan_combo in glycan_combos_list:
        if combo_count == 4:
            glycan_combo_name = intact_chain_combo_label + '|+ ' + \
                                glycan_combo[0] + ' + ' + \
                                glycan_combo[1] + ' + ' + \
                                glycan_combo[2] + ' + ' + \
                                glycan_combo[3]
        else:
            glycan_combo_name = intact_chain_combo_label + '|+ ' + \
                                glycan_combo[0] + ' + ' + \
                                glycan_combo[1]

        total_comp_ctr = Counter(intact_comp_ctr)

        for cur_index in range(combo_count):
            total_comp_ctr += glycan_comps_dict[glycan_combo[cur_index]]

        intact_glycan_comp_dict[glycan_combo_name] = total_comp_ctr

    return intact_glycan_comp_dict
    # End of function calc_intact_glycan_comp_dict()


def gen_intact_glycans_comp_dicts(params_dict: dict,
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
                {dict} -- Intact compositions dictionary.
                {dict} -- Intact compositions with glycosylation dictionary.
                {list} -- Combinations of heavy chains and light chanins.
    '''
    intact_comps_dict: Dict[str, Counter] = {}
    intact_glycan_comps_dict: Dict[str, Counter] = {}

    hc12_lc12_combos: List = cal.hc_lc_combos(hc_chain_list=unique_hc_list,
                                              lc_chain_list=unique_lc_list)

    for cur_combo in hc12_lc12_combos:
        # intact_comps_dict.
        intact_combo_label: str = f'{cur_combo[0]}/{cur_combo[1]}:{cur_combo[2]}/{cur_combo[3]}'
        intact_comp_ctr: Counter = calc_intact_comp_ctr(hc_pos_1_full_reduced_comp_ctr=Counter(reduced_comp_dict[cur_combo[0]]['full_reduced_comp_ctr']),
                                                        hc_pos_2_full_reduced_comp_ctr=Counter(reduced_comp_dict[cur_combo[1]]['full_reduced_comp_ctr']),
                                                        lc_pos_1_full_reduced_comp_ctr=Counter(reduced_comp_dict[cur_combo[2]]['full_reduced_comp_ctr']),
                                                        lc_pos_2_full_reduced_comp_ctr=Counter(reduced_comp_dict[cur_combo[3]]['full_reduced_comp_ctr']),
                                                        total_disulfides=int(params_dict['Total_Disulfides']))
        intact_comps_dict[intact_combo_label] = intact_comp_ctr

        # intact_glycan_comps_dict.
        cur_intact_glycan_comp_dict: Dict[str, Counter] = calc_intact_glycan_comp_dict(intact_comp_ctr=intact_comp_ctr,
                                                                                       intact_chain_combo_label=intact_combo_label,
                                                                                       is_lc_glycos=bool(params_dict['LC_Glycos']),
                                                                                       glycan_comps_dict=hc_lc_glycans_comp_dict)
        intact_glycan_comps_dict.update(cur_intact_glycan_comp_dict)

    return (intact_comps_dict, intact_glycan_comps_dict, hc12_lc12_combos)
    # End of function gen_intact_glycans_comp_dicts
