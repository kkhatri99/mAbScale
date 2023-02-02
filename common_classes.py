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
    Common classes.
'''
# XXX: Marker to show Pylint & Pylance are done checking the code.


# Built-in modules.
from collections import Counter


class ReducedComp:
    '''
        Container for reduced composition results.
    '''

    def __init__(self):
        self.part_reduced_comp_ctr: Counter
        self.part_reduced_msg: str = ''
        self.full_reduced_comp_ctr: Counter
    # End of class ReducedComp


class Chain():
    '''
        Standard values for chains.
    '''
    HC: str = 'HC'
    LC: str = 'LC'
    # End of class Chain()
