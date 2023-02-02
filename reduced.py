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
    Functions reduced compositions.
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built-in modules.
from collections import Counter
from typing import Dict, List

# App modules.
import calc as cal
import common_classes as cc


def calc_hc_reduced_comp(hc_sequence: str,
                         hc_comp_ctr: Counter,
                         cyclize: bool,
                         lys_clip: bool,
                         hc_disulfides: int) -> cc.ReducedComp:
    '''
        Calculate the partial and full reduced compositions for the heavy chain.

        Arguements:
            hc_sequence {str} -- Heavy chain sequence.

            hc_comp_ctr {Counter} -- Heavy chain composition.

            cyclize {bool} -- Is heavy chain cyclized.

            lys_clip {bool} -- Is heavy chain lys clipped.

            hc_disulfides {int} -- Total number of heavy chain disulfides

        Returns:
            {ReducedComps} -- Partial and full reduced compositions for heavy chain.

        Note:
            All cysteines are reduced to begin with (contain a hydrogen/free thiol) these
            will be removed eventually for those in diS bonds.
    '''
    molecule_comp_dict: Dict[str, Counter] = cal.gen_molecule_comp_dict()
    nitrogen: Counter = molecule_comp_dict['Nitrogen']
    hydrogen: Counter = molecule_comp_dict['Hydrogen']
    water: Counter = molecule_comp_dict['Water']
    lys: Counter = molecule_comp_dict['Lys']

    cycliz_msg: str = 'None'
    hc_mod_comp_ctr: Counter = hc_comp_ctr

    # Handle N-Term Cyclization / Pyro_Q or Pyro_E.
    if cyclize:
        if hc_sequence[0] == 'Q':
            hc_mod_comp_ctr -= (cal.mult(1, nitrogen) + cal.mult(3, hydrogen))
            cycliz_msg = 'Pyro Q'
        elif hc_sequence[0] == 'E':
            hc_mod_comp_ctr -= water
            cycliz_msg = 'Pyro E'

    lys_clip_msg: str = 'No'

    # Handle Lysine Clipping.
    if lys_clip and hc_sequence.endswith('K'):
        hc_mod_comp_ctr -= lys
        lys_clip_msg = 'Yes'

    # Remove disulfides.
    # Fully reduced heavy chain, no inter- or intra-chain disulfides.
    hc_full_reduced_comp_ctr = hc_mod_comp_ctr
    # Partially reduced heavy chain, only intra-chain disulfides present
    hc_part_reduced_comp_ctr = hc_mod_comp_ctr - cal.mult(hc_disulfides, cal.mult(2, hydrogen))

    summary_message: str = f'N-Terminal Cyclization: {cycliz_msg} | C-Terminal Lysine Clipping: {lys_clip_msg}'

    reduced_comp = cc.ReducedComp()
    reduced_comp.part_reduced_comp_ctr = hc_part_reduced_comp_ctr
    reduced_comp.part_reduced_msg = summary_message
    reduced_comp.full_reduced_comp_ctr = hc_full_reduced_comp_ctr

    return reduced_comp
    # End of function calc_hc_reduced_comp()


def calc_lc_reduced_comp(lc_comp_ctr: Counter,
                         lc_disulfides: int) -> cc.ReducedComp:
    '''
        Calculate partial and full reduced compositions for the light chain.

        Arguments:
            lc_comp_ctr {Counter} -- Light chain composition.

            lc_disulfides {int} -- Count of light chain disulfides.

        Returns:
            {ReducedComps} -- Various light chain reduced compositions and message.
    '''
    molecule_comp_dict: Dict[str, Counter] = cal.gen_molecule_comp_dict()
    hydrogen: Counter = molecule_comp_dict['Hydrogen']

    lc_full_reduced_comp_ctr: Counter = lc_comp_ctr    # fully reduced light chain, no inter- or intra-chain disulfides
    lc_part_reduced_comp_ctr: Counter = lc_comp_ctr - (cal.mult(lc_disulfides, cal.mult(2, hydrogen)))    # partially-reduced LC, only intra-chain disulfides present

    reduced_comp = cc.ReducedComp()
    reduced_comp.part_reduced_comp_ctr = lc_part_reduced_comp_ctr
    reduced_comp.part_reduced_msg = ''
    reduced_comp.full_reduced_comp_ctr = lc_full_reduced_comp_ctr

    return reduced_comp
    # End of function calc_lc_reduced_comp()


def calc_reduced_glycan_comp_dict(chain_prefix: str,
                                  chain_part_reduced_comp_ctr: Counter,
                                  glycans_comp_dict: Dict[str, Counter]) -> Dict[str, Counter]:
    '''
        Generate a dictionary of reduced chain (heavy or light) with glycan compositions.

        Arguments:
            chain_prefix {str} -- Prefix for chain
            ex: 'HC'

            chain_part_reduced_comp_ctr {Counter} -- Composition of partically reduced chain.

            glycans_comp_dict {Dict[str, Counter]} -- Dictionary of glycans compositions.

        Returns:
            {Dict} -- Dictionary of reduced chain with glycan compositions.
    '''
    chain_prefix += '_'
    glycan_reduced_dict: Dict[str, Counter] = {}

    for glycan_name, glycan_ctr in glycans_comp_dict.items():
        glycan_reduced_dict[chain_prefix + glycan_name] = chain_part_reduced_comp_ctr + \
                                                          Counter(glycan_ctr)

    return glycan_reduced_dict
    # End of function calc_reduced_glycan_comp_dict()


def gen_reduced_comps_dict(chains_dict: dict,
                           params_dict: dict) -> dict:
    '''
        Creates a dictionary of reduced compositions properites.

        Arguments:
            chains_dict {dict} -- Chain properties dictionary.

            params_dict {dict} -- Parameters dictionary.

        Returns:
            {dict} -- Reduced compositions properites dictionary.
            ex: {
                 'HC-1': {'part_reduced_comp_ctr': Counter({'Carbon': 15, 'Hydrogen': 15, 'Oxygen': 6, 'Nitrogen': 3}),
                          'part_reduced_msg': 'N-Terminal Cyclization: None | C-Terminal Lysine Clipping: No',
                          'full_reduced_comp_ctr': Counter({'Hydrogen': 23, 'Carbon': 15, 'Oxygen': 6, 'Nitrogen': 3})},
                 'HC-2': ...
                }
    '''
    reduced_comp_dict: Dict = {}
    chain_num_list: List[str] = ['1', '2']

    # HC Reduced Compositions.
    for chain_num in chain_num_list:
        comp_results: cc.ReducedComp = calc_hc_reduced_comp(hc_sequence=chains_dict[f'HC-{chain_num}']['Sequence'],
                                                            hc_comp_ctr=chains_dict[f'HC-{chain_num}']['Comp'],
                                                            cyclize=bool(params_dict[f'Cyclize_{chain_num}']),
                                                            lys_clip=bool(params_dict[f'Lys_Clip_{chain_num}']),
                                                            hc_disulfides=int(params_dict['HC_Disulfides']))
        reduced_comp_dict[f'HC-{chain_num}'] = {'part_reduced_comp_ctr': comp_results.part_reduced_comp_ctr,
                                                'part_reduced_msg': comp_results.part_reduced_msg,
                                                'full_reduced_comp_ctr': comp_results.full_reduced_comp_ctr}

    # LC Reduced Compositions.
    for chain_num in chain_num_list:
        comp_results: cc.ReducedComp = calc_lc_reduced_comp(lc_comp_ctr=chains_dict[f'LC-{chain_num}']['Comp'],
                                                            lc_disulfides=int(params_dict['LC_Disulfides']))
        reduced_comp_dict[f'LC-{chain_num}'] = {'part_reduced_comp_ctr': comp_results.part_reduced_comp_ctr,
                                                'part_reduced_msg': '',
                                                'full_reduced_comp_ctr': comp_results.full_reduced_comp_ctr}

    return reduced_comp_dict
    # End of fucntion gen_reduced_comps_dict()


def gen_reduced_glycan_comps_dict(chains_dict: dict,
                                  params_dict: dict,
                                  reduced_comp_dict: dict,
                                  hc_glycans_comp_dict: dict,
                                  lc_glycans_comp_dict) -> dict:
    '''
        Creates dictionary of reduced with glycosylation dictionary.

        Arguments:
            chains_dict {dict} -- Chain properties dictionary.

            params_dict {dict} -- Parameters dictionary.

            reduced_comp_dict {dict} -- Reduced compositions dictionary.

            hc_glycans_comp_dict {dict} -- Heavy chain glycan compositions dictionary.

            lc_glycans_comp_dict {dict} -- Light chain glycan compositions dictionary.

        Returns:
            {dict} -- Reduced with glycosylation dictionary.
            ex: {
                'HC-1_Chitobiose core': Counter({'Hydrogen': 41, 'Carbon': 31, 'Oxygen': 16, 'Nitrogen': 5}),
                'HC-1_G0F': Counter({'Hydrogen': 107, 'Carbon': 71, 'Oxygen': 45, 'Nitrogen': 7}),
                'HC-1_G1F': Counter({'Hydrogen': 117, 'Carbon': 77, 'Oxygen': 50, 'Nitrogen': 7}),
                ...
                }
    '''
    reduced_glycos_dict: dict = {}

    # Determine which chains to use.
    reduced_glycos_chains_list: List[str] = ['HC-1']

    if chains_dict['HC-2']['Sequence'] != chains_dict['HC-1']['Sequence']:
        reduced_glycos_chains_list.append('HC-2')

    if bool(params_dict['LC_Glycos']):
        reduced_glycos_chains_list.append('LC-1')

        if chains_dict['LC-2']['Sequence'] != chains_dict['LC-1']['Sequence']:
            reduced_glycos_chains_list.append('LC-2')

    for chain in reduced_glycos_chains_list:
        if chain.startswith('H'):
            glycans_dict = hc_glycans_comp_dict
        else:
            glycans_dict = lc_glycans_comp_dict

        reduced_glycos_dict.update(calc_reduced_glycan_comp_dict(chain_prefix=chain,
                                                                 chain_part_reduced_comp_ctr=reduced_comp_dict[chain]['part_reduced_comp_ctr'],
                                                                 glycans_comp_dict=glycans_dict))

    return reduced_glycos_dict
    # End of function gen_reduced_glycan_comps_dict()
