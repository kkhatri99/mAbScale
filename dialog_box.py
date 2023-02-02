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
    Functions and class for dialog box.
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built in modules.
from collections import Counter
import os
import sys
from typing import Dict, List

# Third party modules.
import PySimpleGUI as sg

# App modules.
import common_functions as cf


def are_chars_in_list_valid(check_string: str,
                            valid_list: List[str]) -> bool:
    '''
        Checks if all characters in string are in the list.

        Arguments:
            check_string {str} -- String with characters to check.

            valid_list {list} -- List of valid values.

        Returns:
            {bool} -- Is string valid.
    '''
    # Convert input parameters to upper case.
    check_string = check_string.strip().upper()
    valid_list = [each_char.upper() for each_char in valid_list]

    valid: bool = True

    for char in check_string:
        if char not in valid_list:
            valid = False
            break

    return valid
    # End of function are_chars_in_list_valid()


def build_dialog_box_layout() -> list:
    '''
        Setup the layout for dialog window.

        Arguments:
            None

        Returns:
            {list} -- Layout for dialog window.
    '''
    label_width: int = 20
    str_input_width: int = 50
    int_input_width: int = 10

    layout = [
        [sg.Text(' ' * 82),
         sg.Text('1'),
         sg.Text(' ' * 88),
         sg.Text('2')],
        [sg.Text(text='Heavy Chain Sequence',
                 size=(label_width, 1)),
         sg.Multiline(size=(str_input_width, 14),
                      background_color = 'LightBlue',
                      key='-HC-1_SEQUENCE-'),
         sg.Multiline(size=(str_input_width, 14),
                      background_color = 'light yellow',
                      key='-HC-2_SEQUENCE-')],
        [sg.Text(text='Heavy Chain Chem Mod',
                 size=(label_width, 1)),
                 sg.InputText(size=(str_input_width, 1),
                              background_color = 'LightBlue',
                              key='-HC-1_CHEM_MOD-'),
                 sg.Text(' '),
                 sg.InputText(size=(str_input_width, 1),
                              background_color = 'light yellow',
                              key='-HC-2_CHEM_MOD-')],
        [sg.Text(text='Light Chain Sequence',
                 size=(label_width, 1)),
         sg.Multiline(size=(str_input_width, 8),
                            background_color = 'SkyBlue1',
                            key='-LC-1_SEQUENCE-'),
         sg.Multiline(size=(str_input_width, 8),
                            background_color = 'LightGoldenrod1',
                            key='-LC-2_SEQUENCE-')],
        [sg.Text(text='Light Chain Chem Mod',
                 size=(label_width, 1)),
                 sg.InputText(size=(str_input_width, 1),
                              background_color = 'SkyBlue1',
                              key='-LC-1_CHEM_MOD-'),
                 sg.Text(' '),
                 sg.InputText(size=(str_input_width, 1),
                              background_color = 'LightGoldenrod1',
                              key='-LC-2_CHEM_MOD-')],

        [sg.Text(' ' * 62),
         sg.Checkbox(text='N-Terminal Cyclization',
                     size=(20, 1),
                     default=True,
                     key='-CYCLIZE_1-'),
         sg.Text(' ' * 43),
         sg.Checkbox(text='N-Terminal Cyclization',
                     size=(20, 1),
                     default=True,
                     key='-CYCLIZE_2-')],
        [sg.Text(' ' * 62),
         sg.Checkbox(text='C-Terminal Lys Clipping',
                     size=(20, 1),
                     default=True,
                     key='-LYS_CLIP_1-'),
         sg.Text(' ' * 43),
         sg.Checkbox(text='C-Terminal Lys Clipping',
                     size=(20, 1),
                     default=True,
                     key='-LYS_CLIP_2-')],

        [sg.Text(' ')],

        [sg.Text(text='Total Number of Disulfides',
                 size=(label_width, 1)),
         sg.InputText('16',
                      size=(int_input_width, 1),
                      justification='right',
                      key='-TOTAL_DISULFIDES-')],      # Used to calculate the intact mass.

        [sg.Text(text=' ')],

        [sg.Text(text='Reduced Analysis')],
        [sg.Text(text='Unreduced HC Disulfides',
                 size=(label_width, 1)),
         sg.InputText('4',
                      size=(int_input_width, 1),
                      justification='right',
                      key='-HC_DISULFIDES-')],         # Disulfides left after reduction in HC.
        [sg.Text(text='Unreduced LC Disulfides',
                 size=(label_width, 1)),
         sg.InputText('2',
                      size=(int_input_width, 1),
                      justification='right',
                      key='-LC_DISULFIDES-')],         # Disulfides left after reduction in LC.

        [sg.Text(' ')],

        [sg.Checkbox(text='Light Chain is Glycosylated',
                     size=(20, 1),
                     default=False,
                     key='-LC_GLYCOS-')],

        [sg.Text(' ')],

        [sg.Text(text='Output Folder',
                 size=(label_width, 1)),
         sg.InputText(size=(95, 1),
                      key='-OUTPUT_FOLDER-'),
         sg.FolderBrowse(target='-OUTPUT_FOLDER-')],
        [sg.Text(text='Excel File Name (no ext)',
                 size=(label_width, 1)),
         sg.InputText(size=(str_input_width, 1),
                      key='-FILE_NAME-'),
         sg.Text(text='.xlsx')],

        [sg.Text(' ')],

        [sg.Text(' ' * 95),
         sg.Submit(),
         sg.Text(' ' * 4),
         sg.Cancel()]
    ]

    return layout
    # End of funciton build_dialog_box_layout()


class ErrorResult:
    '''
        Container for error checking.
    '''

    def __init__(self):
        self.is_error: bool = False
        self.error_list: list = []
    # End of class ErrorResult



def is_int(value: str) -> bool:
    '''
        Check if given string is a integer.

        Arguments:
            value {str} -- String to check.

        Returns:
            {bool} -- If value is an integer.
    '''
    try:
        int(value)
        return True
    except ValueError:
        return False
    # End of function is_int()


def is_in_range_int(num: int,
                    min_val: int,
                    max_val: int) -> bool:
    '''
        Checks if given number is in given range.

        Arguments:
            num {str} -- Value to check.
            min_val {int} -- Minimum valid value.
            max_val {int} -- Maximum valid value.

        Returns:
            {bool} -- Number within range.
    '''
    return min_val <= num <= max_val
    # End of function is_in_range_int()


def is_even(num: int) -> bool:
    '''
        Checks if given number is even.

        Arguments:
            num {int} -- Number to check.

        Returns:
            {bool} -- If number is even.
    '''
    return (num % 2) == 0
    # End of function is_even()


def check_input_params(params_dict: Dict) -> ErrorResult:
    '''
        Checks if the input parameters are valid based on rules.

        Arguments:
            params_dict {Dict} -- Dictionary of input parameters.

        Returns:
            {ErrorResult} -- Error results.
    '''
    aa1_aa3_dict: Dict[str, str] = cf.gen_amino_acid_dict()
    key_list: list = list(aa1_aa3_dict.keys())

    error_list: list = []


    # 1. Check HC-1 and LC-1 have values.
    chain_1_list: list = [{'param': 'HC-1_Sequence', 'label': 'Heavy Chain 1 Sequence'},
                          {'param': 'LC-1_Sequence', 'label': 'Light Chain 1 Sequence'}]

    for chain in chain_1_list:
        if not params_dict[chain['param']]:
            error_list.append(f'"{chain["label"]}" is required.')


    # 2. Check HC-1, LC-1, HC-2, LC-2 sequences contains valid amino acids characters.
    chain_list: list = [{'param': 'HC-1_Sequence', 'label': 'Heavy Chain 1 Sequence'},
                        {'param': 'LC-1_Sequence', 'label': 'Light Chain 1 Sequence'},
                        {'param': 'HC-2_Sequence', 'label': 'Heavy Chain 2 Sequence'},
                        {'param': 'LC-2_Sequence', 'label': 'Light Chain 2 Sequence'}]

    for chain in chain_list:
        if not are_chars_in_list_valid(check_string=params_dict[chain['param']],
                                       valid_list=key_list):
            error_list.append(f'"{chain["label"]}" contains invalid amino acid characters.')


    # 3. Check HC-1, LC-1, HC-2, LC-2 sequences contains no more than 100 amino acids.
    aa_limit: int = 1_000

    for chain in chain_list:
        if len(params_dict[chain['param']]) > aa_limit:
            error_list.append(f'"{chain["label"]}" length is over {aa_limit} amino acids.')


    # 4. Check Chemical Mofidication for HC-1, HC-2, LC-1, LC-2 contain valid elements symbols.
    chain_list: list = [{'param': 'HC-1_Chem_Mod', 'label': 'Heavy Chain 1 Chemical Modificaiton'},
                        {'param': 'LC-1_Chem_Mod', 'label': 'Light Chain 1 Chemical Modificaiton'},
                        {'param': 'HC-2_Chem_Mod', 'label': 'Heavy Chain 2 Chemical Modificaiton'},
                        {'param': 'LC-2_Chem_Mod', 'label': 'Light Chain 2 Chemical Modificaiton'}]

    ele_symbol_mass_dict: Dict = cf.gen_ele_symbol_mass_dict()
    all_ele_sym_list: List[str] = list(ele_symbol_mass_dict.keys())
    all_ele_sym_set: set = set(all_ele_sym_list)

    for chain in chain_list:
        if params_dict[chain['param']]:
            cur_chem_mod_ctr: Counter = cf.formula_to_composition(mol_formula=params_dict[chain['param']])
            param_ele_sym_list: list = list(cur_chem_mod_ctr.keys())
            invalid_ele_sym_list = list(set(param_ele_sym_list) - all_ele_sym_set)
            invalid_ele_sym_list.sort()

            if invalid_ele_sym_list:
                errors_out: str = ','.join(invalid_ele_sym_list)
                error_list.append(f'"{chain["label"]}" contains invalid invalid element sysmbol(s): {errors_out}.')


    # 5. Check Total Disulfides, HC Disulfides, and LC Disulfides are:
    #      a. Integers
    #      b. Within valid range
    #      c. Even
    min_value_check: int = 0
    max_value_check: int = 100

    disulfides_list: list = [{'param': 'Total_Disulfides', 'label': 'Total Number of Disulfides'},
                             {'param': 'HC_Disulfides', 'label': 'Unreduced HC Disulfides'},
                             {'param': 'LC_Disulfides', 'label': 'Unreduced LC Disulfides'}]

    for disulfide in disulfides_list:
        if not is_int(params_dict[disulfide['param']]):
            error_list.append(f'"{disulfide["label"]}" is not an integer.')
        else:
            if not is_in_range_int(num=int(params_dict[disulfide['param']]),
                                   min_val=min_value_check,
                                   max_val=max_value_check):
                error_list.append(f'"{disulfide["label"]}" is not within valid range ({min_value_check} to {max_value_check}).')
            else:
                if not is_even(num=int(params_dict[disulfide['param']])):
                    error_list.append(f'"{disulfide["label"]}" is not an even number.')


    # 6. Check Output Folder exists.
    is_dir = os.path.isdir(params_dict['Output_Folder'])
    if not is_dir:
        error_list.append('"Output Folder" is not valid')


    # 7. Check Folder Name has a value.
    if not params_dict['File_Name']:
        error_list.append('"File Name" is required')


    current_error = ErrorResult()

    if error_list:
        current_error.is_error = True
        current_error.error_list = error_list

    return current_error
    # End of function check_input_params()


def show_dialog_box() -> Dict:
    '''
        Show dialog box for parameters.

        Arguments:
            None

        Returns:
            {Dict} -- Dictionary of parameters.
    '''
    params_dict: Dict = {}

    sg.change_look_and_feel('Reddit')
    layout = build_dialog_box_layout()

    window = sg.Window(title='mAbScale',
                       icon='mass_calc.ico',
                       layout=layout)
    values: dict = {}

    while True:
        event, values = window.read()       # type: ignore

        if event in (sg.WIN_CLOSED, 'Cancel'):
            sys.exit("0")

        params_dict['HC-1_Sequence'] = values['-HC-1_SEQUENCE-'].strip().upper()
        params_dict['HC-1_Chem_Mod'] = values['-HC-1_CHEM_MOD-'].strip()

        params_dict['LC-1_Sequence'] = values['-LC-1_SEQUENCE-'].strip().upper()
        params_dict['LC-1_Chem_Mod'] = values['-LC-1_CHEM_MOD-'].strip()

        params_dict['Cyclize_1'] = values['-CYCLIZE_1-']
        params_dict['Lys_Clip_1'] = values['-LYS_CLIP_1-']

        params_dict['HC-2_Sequence'] = values['-HC-2_SEQUENCE-'].strip().upper()
        params_dict['HC-2_Chem_Mod'] = values['-HC-2_CHEM_MOD-'].strip()

        params_dict['LC-2_Sequence'] = values['-LC-2_SEQUENCE-'].strip().upper()
        params_dict['LC-2_Chem_Mod'] = values['-LC-2_CHEM_MOD-'].strip()

        params_dict['Cyclize_2'] = values['-CYCLIZE_2-']
        params_dict['Lys_Clip_2'] = values['-LYS_CLIP_2-']

        # Account for mAb shorthand entry (Entry of just HC-1 and LC-1 respresents a mAb).
        if not params_dict['HC-2_Sequence'] and not params_dict['LC-2_Sequence']:
            params_dict['HC-2_Sequence'] = params_dict['HC-1_Sequence']
            params_dict['HC-2_Chem_Mod'] = params_dict['HC-1_Chem_Mod']

            params_dict['LC-2_Sequence'] = params_dict['LC-1_Sequence']
            params_dict['LC-2_Chem_Mod'] = params_dict['LC-1_Chem_Mod']

            params_dict['Cyclize_2'] = params_dict['Cyclize_1']
            params_dict['Lys_Clip_2'] = params_dict['Lys_Clip_1']
        else:   # Bispecific.
            if not params_dict['HC-2_Sequence']:
                params_dict['HC-2_Sequence'] = params_dict['HC-1_Sequence']

            if not params_dict['LC-2_Sequence']:
                params_dict['LC-2_Sequence'] = params_dict['LC-1_Sequence']

        params_dict['Total_Disulfides'] = values['-TOTAL_DISULFIDES-'].strip() or 0

        params_dict['HC_Disulfides'] = values['-HC_DISULFIDES-'].strip() or 0
        params_dict['LC_Disulfides'] = values['-LC_DISULFIDES-'].strip() or 0

        params_dict['LC_Glycos'] = values['-LC_GLYCOS-']

        params_dict['Output_Folder'] = values['-OUTPUT_FOLDER-'].strip()
        params_dict['File_Name'] = values['-FILE_NAME-'].strip()
        params_dict['File_Path'] = params_dict['Output_Folder'] + '/' + params_dict['File_Name'] + '.xlsx'

        current_error = check_input_params(params_dict=params_dict)

        if not current_error.is_error:
            break

        error_count: int = len(current_error.error_list)
        error_list_counter: list = [' ']
            # Add counter to each error.
        for counter, err_msg in enumerate(current_error.error_list, start=1):
            error_list_counter.append(f'{counter:2d}) {err_msg}')

        error_list_counter.append(' ')
        error_msg: str = "\n".join(error_list_counter)

        sg.popup_ok(error_msg,
                    icon='mass_calc.ico',
                    title=f'{error_count} Error(s) Found',
                    background_color='pale violet red',
                    line_width=230)
    window.close()

    return params_dict
    # End of function show_dialog_box()
