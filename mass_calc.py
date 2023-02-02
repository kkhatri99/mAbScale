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
    mAbScale Mass Calculator.

    Authors: Kshitij Khatri
             Thomas Harkins
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built-in modules.
from collections import Counter
from datetime import timedelta
from timeit import default_timer as timer
from typing import Dict, List

# Third party modules.
import PySimpleGUI as sg

# App modules.
import calc as cal
import common_classes as cc
import dialog_box as db
import half_body as hb
import intact as i
import reduced as red
import report as rep


def main() -> None:
    '''
        Main function.

        Arguments:
            None

        Returns:
            Nothing
    '''
    start = timer()
    print('Start of mass calculations')


    params_dict: Dict = db.show_dialog_box()

    hclc_chain_list: List[str] = ['HC-1', 'LC-1', 'HC-2', 'LC-2']
    unique_hc_chain_list: List[str] = cal.gen_unique_chain_list(chain_type=cc.Chain.HC,
                                                                params_dict=params_dict)
    unique_lc_chain_list: List[str] = cal.gen_unique_chain_list(chain_type=cc.Chain.LC,
                                                                params_dict=params_dict)
    chains_dict: dict = cal.build_chain_dict(chains_list=hclc_chain_list,
                                             params_dict=params_dict)

    hc_glycans_comp_dict: Dict[str, Counter] = cal.gen_glycan_comp_dict(chain_type=cc.Chain.HC)
    lc_glycans_comp_dict: Dict[str, Counter] = cal.gen_glycan_comp_dict(chain_type=cc.Chain.LC)
    hc_lc_glycans_comp_dict: Dict[str, Counter] = {**hc_glycans_comp_dict, **lc_glycans_comp_dict}

    chem_mod_dict: dict = cal.build_chem_mod_dict(chains_list=hclc_chain_list,
                                                  params_dict=params_dict)

    reduced_comps_dict: dict = red.gen_reduced_comps_dict(chains_dict=chains_dict,
                                                          params_dict=params_dict)
    reduced_glycan_comps_dict: dict = red.gen_reduced_glycan_comps_dict(chains_dict=chains_dict,
                                                                        params_dict=params_dict,
                                                                        reduced_comp_dict=reduced_comps_dict,
                                                                        hc_glycans_comp_dict=hc_glycans_comp_dict,
                                                                        lc_glycans_comp_dict=lc_glycans_comp_dict)

    intact_comps_dict, intact_glycan_comps_dict, intact_combo_list = i.gen_intact_glycans_comp_dicts(params_dict=params_dict,
                                                                                                     reduced_comp_dict=reduced_comps_dict,
                                                                                                     unique_hc_list=unique_hc_chain_list,
                                                                                                     unique_lc_list=unique_lc_chain_list,
                                                                                                     hc_lc_glycans_comp_dict=hc_lc_glycans_comp_dict)
    intact_chem_mod_dict: dict = cal.build_chem_mod_chain_combo_dict(chain_combo_list=intact_combo_list,
                                                                     num_of_chains=4,
                                                                     chem_mod_dict=chem_mod_dict)

    halfbody_comps_dict: dict = {}
    halfbody_glycan_comps_dict: Dict[str, Counter] = {}
    halfbody_chem_mod_dict: Dict[str, Counter] = {}

    is_hetrodimer: bool = len(unique_hc_chain_list) != 1 or len(unique_lc_chain_list) != 1

    if is_hetrodimer:
        halfbody_comps_dict, halfbody_glycan_comps_dict, halfbody_combo_list = hb.gen_halfbody_glycans_comp_dicts(params_dict=params_dict,
                                                                                                                  reduced_comp_dict=reduced_comps_dict,
                                                                                                                  unique_hc_list=unique_hc_chain_list,
                                                                                                                  unique_lc_list=unique_lc_chain_list,
                                                                                                                  hc_lc_glycans_comp_dict=hc_glycans_comp_dict)
        halfbody_chem_mod_dict = cal.build_chem_mod_chain_combo_dict(chain_combo_list=halfbody_combo_list,
                                                                           num_of_chains=2,
                                                                           chem_mod_dict=chem_mod_dict)

    rep.write_results_to_file(params_dict=params_dict,
                              hc_glycans_comp_dict=hc_glycans_comp_dict,
                              lc_glycans_comp_dict=lc_glycans_comp_dict,
                              is_hetrodimer=is_hetrodimer,
                              chem_mod_dict=chem_mod_dict,
                              intact_chem_mod_dict=intact_chem_mod_dict,
                              halfbody_chem_mod_dict=halfbody_chem_mod_dict,
                              reduced_comps_dict=reduced_comps_dict,
                              intact_comps_dict=intact_comps_dict,
                              halfbody_comps_dict=halfbody_comps_dict,
                              reduced_glycan_comps_dict=reduced_glycan_comps_dict,
                              intact_glycan_comps_dict=intact_glycan_comps_dict,
                              halfbody_glycan_comps_dict=halfbody_glycan_comps_dict)


    print()
    print('End of mass calculations')

    end = timer()
    elapsed_time = timedelta(seconds=end-start)
    print()
    print(f'Elapsed time: {elapsed_time}')

    sg.PopupAutoClose('Mass Calc\nReport created',
                      background_color='green',
                      text_color='white',
                      auto_close_duration=3,
                      button_type=5,
                      no_titlebar=True,
                      font='Helvetica 14')
    # End of function main()


if __name__ == '__main__':
    main()
