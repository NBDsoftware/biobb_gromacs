#!/usr/bin/env python3

"""Module containing the Ndx2resttop class and the command line interface."""
import fnmatch
import argparse
from typing import Optional
from typing import Any
from pathlib import Path
from biobb_common.generic.biobb_object import BiobbObject
from biobb_common.configuration import settings
from biobb_common.tools import file_utils as fu
from biobb_common.tools.file_utils import launchlogger


class Ndx2resttop(BiobbObject):
    """
    | biobb_gromacs Ndx2resttop
    | Generate a restrained topology from an index NDX file.
    | This module automatizes the process of restrained topology generation starting from an index NDX file.

    Args:
        input_ndx_path (str): Path to the input NDX index file. File type: input. `Sample file <https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs_extra/ndx2resttop.ndx>`_. Accepted formats: ndx (edam:format_2033).
        input_top_zip_path (str): Path the input TOP topology in zip format. File type: input. `Sample file <https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/data/gromacs_extra/ndx2resttop.zip>`_. Accepted formats: zip (edam:format_3987).
        output_top_zip_path (str): Path the output TOP topology in zip format. File type: output. `Sample file <https://github.com/bioexcel/biobb_gromacs/raw/master/biobb_gromacs/test/reference/gromacs_extra/ref_ndx2resttop.zip>`_. Accepted formats: zip (edam:format_3987).
        properties (dict - Python dictionary object containing the tool parameters, not input/output files):
            * **posres_names** (*str*) - ("CUSTOM_POSRES") Space-separated preprocessor macro names, one per triplet (e.g. ``"CHAIN_A_POSRES CHAIN_B_POSRES"``). Activated at grompp time via ``-D <name>``. Must match the length of ref_rest_mol_triplet_list.
            * **force_constants** (*str*) - ("500 500 500") Three space-separated force constants Fx Fy Fz in kJ/mol/nm². Written verbatim into each ITP line so the same value is applied to every restrained atom.
            * **ref_rest_mol_triplet_list** (*str*) - (None) Comma-separated triplets ``(reference_group, restrain_group, molecule_name)``, one per molecule to restrain.``reference_group``: NDX group containing *all* atoms of the molecule in global (system-wide) GROMACS numbering. Used only as a coordinate frame to convert global indices to molecule-local 1-based indices. Example: ``r_1-9280``. ``restrain_group``: NDX group with the subset of atoms to actually restrain (must be a strict subset of reference_group). Example: ``P_&_r_1-9280``. ``molecule_name``: the ``[ moleculetype ]`` name in the topology used to locate where to splice the ``#ifdef`` block. Example: ``1TSR``.

    Examples:
        This is a use example of how to use the building block from Python::

            from biobb_gromacs.gromacs_extra.ndx2resttop import ndx2resttop
            prop = { 'ref_rest_mol_triplet_list': '( Chain_A, Chain_A_noMut, Protein_chain_A ), ( Chain_B, Chain_B_noMut, Protein_chain_B ), ( Chain_C, Chain_C_noMut, Protein_chain_C ), ( Chain_D, Chain_D_noMut, Protein_chain_D )' }
            ndx2resttop(input_ndx_path='/path/to/myIndex.ndx',
                        input_top_zip_path='/path/to/myTopology.zip',
                        output_top_zip_path='/path/to/newTopology.zip',
                        properties=prop)

    Info:
        * wrapped_software:
            * name: In house
            * license: Apache-2.0
        * ontology:
            * name: EDAM
            * schema: http://edamontology.org/EDAM.owl
    """

    def __init__(self, input_ndx_path: str, input_top_zip_path: str, output_top_zip_path: str,
                 properties: Optional[dict] = None, **kwargs) -> None:
        properties = properties or {}

        # Call parent class constructor
        super().__init__(properties)
        self.locals_var_dict = locals().copy()

        # Input/Output files
        self.io_dict = {
            "in": {"input_ndx_path": input_ndx_path, "input_top_zip_path": input_top_zip_path},
            "out": {"output_top_zip_path": output_top_zip_path}
        }

        # Properties specific for BB
        self.posres_names = properties.get('posres_names')
        self.force_constants = properties.get('force_constants', '500 500 500')
        self.ref_rest_mol_triplet_list = properties.get('ref_rest_mol_triplet_list')

        # Check the properties
        self.check_properties(properties)
        self.check_arguments()

    @launchlogger
    def launch(self) -> int:
        """Execute the :class:`Ndx2resttop <gromacs_extra.ndx2resttop.Ndx2resttop>` object."""
        # Setup Biobb
        if self.check_restart():
            return 0

        top_file = fu.unzip_top(zip_file=self.io_dict['in'].get("input_top_zip_path", ""), out_log=self.out_log)

        # Create group dictionary from the ndx index file
        # Each group will have a list with its starting and ending lines in the ndx file
        groups_dic: dict[str, Any] = {}
        current_group = ''
        previous_group = ''
        
        # Iterate over the lines of the ndx file
        ndx_lines = open(self.io_dict['in'].get("input_ndx_path", "")).read().splitlines()
        for index, line in enumerate(ndx_lines):
            
            # Found a new group in the index file
            if line.startswith('['):
                current_group = line
                
                # Set the starting line of the new group
                groups_dic[current_group] = [index, 0]

                # Check if there was a previous group to set its ending line
                if previous_group != '':
                    groups_dic[previous_group][1] = index
                    
                # Update the previous group variable
                previous_group = current_group

            # If we reach the last line of the file, set the ending line of the last group
            if index == len(ndx_lines)-1:
                groups_dic[current_group][1] = index + 1

        # Catch groups with just one line
        for group_name in groups_dic.keys():
            if (groups_dic[group_name][0]+1) == groups_dic[group_name][1]:
                groups_dic[group_name][1] += 1

        fu.log('groups_dic: '+str(groups_dic), self.out_log, self.global_log)

        self.ref_rest_mol_triplet_list = [tuple(elem.strip(' ()').replace(' ', '').split(',')) for elem in str(self.ref_rest_mol_triplet_list).split('),')]
        fu.log('ref_rest_mol_triplet_list: ' + str(self.ref_rest_mol_triplet_list), self.out_log, self.global_log)

        if self.posres_names:
            self.posres_names = [elem.strip() for elem in self.posres_names.split()]
            fu.log('posres_names: ' + str(self.posres_names), self.out_log, self.global_log)
        else:
            self.posres_names = ['CUSTOM_POSRES']*len(self.ref_rest_mol_triplet_list)
            fu.log('posres_names: ' + str(self.posres_names), self.out_log, self.global_log)

        # Check if the number of posres_names matches the number of ref_rest_mol_triplet_list
        if len(self.posres_names) != len(self.ref_rest_mol_triplet_list):
            raise ValueError("If posres_names is provided, it should match the number of ref_rest_mol_triplet_list")

        for triplet, posre_name in zip(self.ref_rest_mol_triplet_list, self.posres_names):

            reference_group, restrain_group, molecule_name = triplet

            fu.log('Reference group: '+reference_group, self.out_log, self.global_log)
            fu.log('Restrain group: '+restrain_group, self.out_log, self.global_log)
            fu.log('Molecule name: '+molecule_name, self.out_log, self.global_log)
            self.io_dict['out']["output_itp_path"] = fu.create_name(path=str(Path(top_file).parent), prefix=self.prefix, step=self.step, name=restrain_group+'_posre.itp')

            # Mapping atoms from absolute enumeration to Chain relative enumeration
            fu.log('reference_group_index: start_closed:'+str(groups_dic['[ '+reference_group+' ]'][0]+1)+' stop_open: '+str(groups_dic['[ '+reference_group+' ]'][1]), self.out_log, self.global_log)
            reference_group_list = [int(elem) for line in ndx_lines[groups_dic['[ '+reference_group+' ]'][0]+1: groups_dic['[ '+reference_group+' ]'][1]] for elem in line.split()]
            fu.log('restrain_group_index: start_closed:'+str(groups_dic['[ '+restrain_group+' ]'][0]+1)+' stop_open: '+str(groups_dic['[ '+restrain_group+' ]'][1]), self.out_log, self.global_log)
            restrain_group_list = [int(elem) for line in ndx_lines[groups_dic['[ '+restrain_group+' ]'][0]+1: groups_dic['[ '+restrain_group+' ]'][1]] for elem in line.split()]
            selected_list = [reference_group_list.index(atom)+1 for atom in restrain_group_list]
            # Creating new ITP with restrains
            with open(self.io_dict['out'].get("output_itp_path", ''), 'w') as f:
                fu.log('Creating: '+str(f)+' and adding the selected atoms force constants', self.out_log, self.global_log)
                f.write('[ position_restraints ]\n')
                f.write('; atom  type      fx      fy      fz\n')
                for atom in selected_list:
                    f.write(str(atom)+'     1  '+self.force_constants+'\n')

            multi_chain = False
            # For multi-chain topologies
            # Including new ITP in the corresponding ITP-chain file
            for file_dir in Path(top_file).parent.iterdir():
                if "posre" not in file_dir.name and not file_dir.name.endswith("_pr.itp"):
                    if fnmatch.fnmatch(str(file_dir), "*"+molecule_name+"*.itp"):
                        multi_chain = True
                        with open(str(file_dir), 'a') as f:
                            fu.log('Opening: '+str(f)+' and adding the ifdef include statement', self.out_log, self.global_log)
                            f.write('\n')
                            f.write('; Include Position restraint file\n')
                            f.write('#ifdef '+str(posre_name)+'\n')
                            f.write('#include "'+str(Path(self.io_dict['out'].get("output_itp_path", "")).name)+'"\n')
                            f.write('#endif\n')
                            f.write('\n')

            # For single-chain topologies
            # Including new ITP in the TOP file
            if not multi_chain:

                # Read all lines of the top file
                with open(top_file, 'r') as f:
                    lines = f.readlines()

                in_moleculetype = False
                found_molecule = False
                index = None

                # Find the line where the custom position restraints are going to be included
                for i, line in enumerate(lines):
                    stripped = line.strip()

                    # Find insertion point: just before the next section or existing restraint block
                    if found_molecule:
                        is_next_section = stripped.startswith('[') and (
                            'system' in stripped or 'molecules' in stripped or 'moleculetype' in stripped)
                        if is_next_section or line.startswith('#include ') or line.startswith('#ifdef POSRES'):
                            index = i
                            break

                    # Detect [ moleculetype ] directive
                    if stripped.startswith('[') and 'moleculetype' in stripped:
                        in_moleculetype = True
                        continue

                    # Find the molecule name line inside the [ moleculetype ] block
                    if in_moleculetype:
                        if stripped and not stripped.startswith(';'):
                            in_moleculetype = False
                            if not stripped.startswith('[') and stripped.split()[0] == molecule_name:
                                found_molecule = True

                if not found_molecule:
                    raise ValueError(f"Molecule type '{molecule_name}' not found in the topology file")

                if index is None:
                    index = len(lines)

                # Include the custom position restraints in the top file
                lines.insert(index, '\n')
                lines.insert(index + 1, '; Include Position restraint file\n')
                lines.insert(index + 2, '#ifdef '+str(posre_name)+'\n')
                lines.insert(index + 3, '#include "'+str(Path(self.io_dict['out'].get("output_itp_path", "")).name)+'"\n')
                lines.insert(index + 4, '#endif\n')
                lines.insert(index + 5, '\n')

                # Write the new top file
                with open(top_file, 'w') as f:
                    f.writelines(lines)

        # zip topology
        fu.zip_top(zip_file=self.io_dict['out'].get("output_top_zip_path", ""), top_file=top_file, out_log=self.out_log, remove_original_files=self.remove_tmp)

        # Remove temporal files
        self.remove_tmp_files()

        self.check_arguments(output_files_created=True, raise_exception=False)
        return 0


def ndx2resttop(input_ndx_path: str, input_top_zip_path: str, output_top_zip_path: str,
                properties: Optional[dict] = None, **kwargs) -> int:
    """Create :class:`Ndx2resttop <gromacs_extra.ndx2resttop.Ndx2resttop>` class and
    execute the :meth:`launch() <gromacs_extra.ndx2resttop.Ndx2resttop.launch>` method."""
    return Ndx2resttop(input_ndx_path=input_ndx_path,
                       input_top_zip_path=input_top_zip_path,
                       output_top_zip_path=output_top_zip_path,
                       properties=properties, **kwargs).launch()


ndx2resttop.__doc__ = Ndx2resttop.__doc__


def main():
    """Command line execution of this building block. Please check the command line documentation."""
    parser = argparse.ArgumentParser(description="Wrapper for the GROMACS extra ndx2resttop module.",
                                     formatter_class=lambda prog: argparse.RawTextHelpFormatter(prog, width=99999))
    parser.add_argument('-c', '--config', required=False, help="This file can be a YAML file, JSON file or JSON string")

    # Specific args of each building block
    required_args = parser.add_argument_group('required arguments')
    required_args.add_argument('--input_ndx_path', required=True)
    required_args.add_argument('--input_top_zip_path', required=True)
    required_args.add_argument('--output_top_zip_path', required=True)

    args = parser.parse_args()
    config = args.config if args.config else None
    properties = settings.ConfReader(config=config).get_prop_dic()

    # Specific call of each building block
    ndx2resttop(input_ndx_path=args.input_ndx_path, input_top_zip_path=args.input_top_zip_path,
                output_top_zip_path=args.output_top_zip_path, properties=properties)


if __name__ == '__main__':
    main()
