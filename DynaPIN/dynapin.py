import argparse
import os
from DynaPIN.Tables import dynapin
from DynaPIN.Plots import Plotter
import json
import time
from datetime import datetime
import warnings
warnings.filterwarnings('ignore')

try:
    from DynaPIN import __version__
except ImportError:
    __version__ = "0.1.5" 

def print_banner():
    """
    Prints a clean, professional banner for the package.
    """
    
    banner_art = r"""
         ____                    ____  _  _   _
        |  _ \ _   _ _ __   __ _|  _ \(_)| \ | |
        | | | | | | | '_ \ / _` | |_) | ||  \| |
        | |_| | |_| | | | | (_| |  __/| || |\  |
        |____/ \__,_|_| |_|\__,_|_|   |_||_| \_|
               _  | | 
               \___/
        """
    print(banner_art)
    
def print_stars(r):
    for i in range(r):
        for j in range(29):
            if (i + j) % 3 == 0:
                print("★", end=' ')
            else:
                print("☆", end=' ')
        print()




parser = argparse.ArgumentParser()

parser.add_argument('-t', '--trajectory_file', type=str, help="Trajectory file in .dcd or .pdb formats. If .dcd, pelase provide topology file with '--topology_file' arguman.")
parser.add_argument('-c', '--commands', type=str, help="Commands to run (comma-separated). Choices:\nTables -> 'all_analysis', 'QualityControl', 'ResidueBased', 'InteractionBased'\nPlots -> 'all_plots', plot-specific names (e.g. 'PlotRMSD') or grouped names 'QualityControl', 'ResidueBased', 'InteractionBased'.")
#job_name
parser.add_argument('-o', '--output_dir', type=str, help='The name of the job, if null, DynaPIN will generate an output directory from input file.')
#tables
parser.add_argument('--foldx_path', type=str, help="(Optional) Path to FoldX executable. If omitted, FoldX calculations are skipped.")
parser.add_argument('--topology_file', type=str, help="Path of topology file for dcd trajectories.")
parser.add_argument('-s', '--stride', type=int, help="Stride value.")
parser.add_argument('-sm', '--split_models', type=bool, help="Whether models will be splitted or not.")
parser.add_argument('-ch', '--chains', type=str, help="For heteromers, you can decide which two chains to be analyzed. Provide input such as 'A,B'.")
parser.add_argument('--rmsd_rs', type=str, help="Path to reference structure for RMSD analyse.") 
parser.add_argument('--rmsd_rf', type=int, help="Reference frame for RMSD analyze.")
parser.add_argument('--table_json', type=str, help="JSON file path for analyse runs.")
parser.add_argument('--plot_json', type=str, help="JSON file path for plot runs.")

args = parser.parse_args()

PLOT_GROUPS = {
    'QualityControl': [
        'PlotRMSD', 'PlotRG', 'PlotRMSF', 'PlotiRMSD',
        'PlotlRMSD', 'PlotDockQ', 'PlotFnonnat', 'PlotFnat'
    ],
    'ResidueBased': [
        'PlotBiophys', 'PlotSASA', 'PlotResEne', 'PlotDSSP'
    ],
    'InteractionBased': [
        'PlotPairwiseFreq'
    ],
}
ALL_PLOT_TASKS = []
for group in PLOT_GROUPS.values():
    for cmd in group:
        if cmd not in ALL_PLOT_TASKS:
            ALL_PLOT_TASKS.append(cmd)
if args.commands:
    commands = args.commands.split(',')


if args.table_json:
    table_commands = []
    a_json_ = True

    f = open(args.table_json)
    table_data = json.load(f)

    if table_data['QualityControl']['Run']:
        table_commands.append('QualityControl')

    if table_data['ResidueBased']['Run']:
        table_commands.append('ResidueBased')

    if table_data['InteractionBased']['Run']:
        table_commands.append('InteractionBased')
    
else:
    a_json_ = False
    if args.commands:
        if 'all_analysis' in commands:
            table_commands = ['QualityControl', 'ResidueBased', 'InteractionBased']
        else:
            table_commands = [x for x in commands if 'plot' not in x.lower()]
    else:
        table_commands = None


def _dedupe_preserve(seq):
    seen = set()
    ordered = []
    for item in seq:
        if item not in seen:
            seen.add(item)
            ordered.append(item)
    return ordered

if args.plot_json:
    plot_commands = []
    p_json = True

    f2 = open(args.plot_json)
    plot_data = json.load(f2)

    per_plot_entries = ['PlotRMSD', 'PlotRG', 'PlotRMSF', 'PlotBiophys', 'PlotSASA',
                        'PlotiRMSD', 'PlotlRMSD', 'PlotDockQ', 'PlotFnonnat',
                        'PlotFnat', 'PlotPairwiseFreq', 'PlotResEne', 'PlotDSSP']

    for entry in per_plot_entries:
        if plot_data.get(entry, {}).get('Run'):
            plot_commands.append(entry)

    for group_name, members in PLOT_GROUPS.items():
        group_data = plot_data.get(group_name)
        if isinstance(group_data, dict) and group_data.get('Run'):
            plot_commands.extend(members)

    plot_commands = _dedupe_preserve(plot_commands)
    if not plot_commands:
        plot_commands = None
    

else:
    p_json = False
    if args.commands:
        plot_commands = []
        if 'all_plots' in commands:
            plot_commands = list(ALL_PLOT_TASKS)
        else:
            for cmd in commands:
                if cmd in PLOT_GROUPS:
                    plot_commands.extend(PLOT_GROUPS[cmd])
                elif 'plot' in cmd.lower():
                    plot_commands.append(cmd)
        plot_commands = _dedupe_preserve(plot_commands)
        if not plot_commands:
            plot_commands = None
    else:
        plot_commands = None
    
def main():
    print("\n")
    print_stars(1)
    print("\n")
    print_banner()
    print("\n")
    print(f" Dynamic Interface Analysis of Protein Complexes (v{__version__})")
    print(f"              Generated by Karaca Lab, IBG")
    print("\n")
    print_stars(1)
    print("\n")
    time.sleep(0.5)
    
    mol = None
    #tables
    if table_commands:

        if a_json_:
            trajectory_file = table_data['trajectory_file']
            output_dir = table_data['output_dir']
            stride = table_data['stride']
            split_models = table_data['split_models']
            chains = table_data['chains']
            topology_file = table_data['topology_file']
        

        else:
            trajectory_file = args.trajectory_file
            output_dir = args.output_dir
            stride = args.stride
            split_models = args.split_models
            chains = args.chains
            topology_file = args.topology_file
        
        if not stride:
            stride = 1
        if not split_models:
            split_models=True
            
        if output_dir:
            output_dir = os.path.abspath(output_dir)
 
        def _paths(p):
            if p is None:
                return None
            if os.path.isabs(p) and os.path.exists(p):
                return p
            if os.path.exists(p):
                return os.path.abspath(p)
            if output_dir:
                cand = os.path.abspath(os.path.join(output_dir, p))
                if os.path.exists(cand):
                    return cand
            return p 
        
        trajectory_file = _paths(trajectory_file)
        topology_file   = _paths(topology_file) if topology_file else None
 
        print('Creating DynaPIN class...\n')
        mol = dynapin(trajectory_file=trajectory_file, stride=stride, split_models=split_models, chains=chains, output_dir=output_dir, topology_file=topology_file
                        )

        print(f'Your DynaPIN Class has been created with the following parameters:\n\tOutput directory:{mol.output_dir}\n\tTrajectory File: {trajectory_file}\n\tTopology File: {topology_file}\n\tStride: {stride}\n\tSplit Models: {split_models}\n\tChain Selection: {chains}\n')
        print_stars(1)
        print("\n")


        if 'QualityControl' in table_commands:

            if a_json_:
                rmsd_rs = table_data['QualityControl']['rmsd_data']['ref_struc']
                rmsd_rf = table_data['QualityControl']['rmsd_data']['ref_frame']

            else:
                rmsd_rs = args.rmsd_rs
                rmsd_rf = args.rmsd_rf

            if not rmsd_rf:
                rmsd_rf = 0
            else:
                rmsd_rf = int(rmsd_rf/stride)
                print(rmsd_rf)


            print('Running Quality Control Analysis...\n')
            start_time = datetime.now()
            mol.run_quality_control(rmsd_data={'ref_struc':rmsd_rs, 'ref_frame':rmsd_rf})
            end_time = datetime.now()
            print(f"Quality Control Analysis has run successfully!\nRunning duration: {end_time - start_time}\n")
            print_stars(1)
            print("\n")
        if 'ResidueBased' in table_commands:

            if a_json_:
                foldx_path = table_data['ResidueBased']['FoldX_path']
                if 'Run_DSSP' in table_data['ResidueBased'] and not table_data['ResidueBased']['Run_DSSP']:
                    warnings.warn("The 'Run_DSSP' option is deprecated; DSSP now always runs during ResidueBased analysis.")
            else:
                foldx_path = args.foldx_path or None
            
            if foldx_path is None:
                warnings.warn("FoldX path not provided. Residue energy calculations will be skipped.")

            print('Running Residue Based Analysis...\n')
            start_time = datetime.now()
            mol.run_res_based(foldx_path=foldx_path)
            end_time = datetime.now()
            print(f"Residue Based Analysis has run successfully!\nRunning duration: {end_time - start_time}\n")
            print_stars(1)
            print("\n")

        if 'InteractionBased' in table_commands:

            if a_json_:
                all_hph = table_data['InteractionBased']['get_all_hph']

            else:
                all_hph = getattr(args, 'all_hph', None)
            
            if all_hph is None or all_hph is False:
                all_hph = False
            
            print('Running Interaction Based Analysis... Running Interaction Analyze may take several minutes to hours, based on your input.\n')
            start_time = datetime.now()
            mol.run_inter_based(get_all_hph=all_hph)
            end_time = datetime.now()
            print(f"Interaction Based Analysis has run successfully!\nRunning duration: {end_time - start_time}\n")
            print_stars(1)
            print("\n")

        if not a_json_:
            mol._get_params_()
            print("Your run parameters are downloaded as table_params.json")
        else:
            f.close()

    #plots
    if plot_commands:
        if p_json:
            output_dir = plot_data['output_dir']
        else:
            # If analysis was just run, use its output dir. Otherwise, use arg.
            output_dir = mol.output_dir if mol else args.output_dir


        print_stars(1)
        print("\n")
        draw = Plotter(output_dir=output_dir)
        print(f'Your Plotter Class has been created\n')

        print_stars(1)
        print("\n")

        print('Your plots are processing...\n')
        print_stars(1)
        print("\n")

        if 'PlotRMSD' in plot_commands:
            
            if p_json:
                rmsd_tpath = plot_data['PlotRMSD']['rmsd_table_path']
                draw.plot_rmsd(path=rmsd_tpath)
            else:
                draw.plot_rmsd()
            print('Time Evolution of backbone RMSD plot is done!\n')

        if 'PlotRG' in plot_commands:
            
            if p_json:
                rg_tpath = plot_data['PlotRG']['rg_table_path']
                draw.plot_rg(path=rg_tpath)
            else:
                draw.plot_rg()
            print('Dynamic Variation of Radius of Gyration (R_g) plot is done!\n')
            
        if 'PlotRMSF' in plot_commands:

            if p_json:
                rmsf_tpath = plot_data['PlotRMSF']['rmsf_table_path']
                rmsf_itpath = plot_data['PlotRMSF']['rmsf_intf_res_table']
                draw.plot_rmsf(rmsf_path=rmsf_tpath, intf_path=rmsf_itpath)
            else:
                draw.plot_rmsf()
            print('Residue-wise backbone RMSF plot is done!\n')
        
        if 'PlotBiophys' in plot_commands:

            if p_json:
                biophys_tpath = plot_data['PlotBiophys']['biophys_table_path']
                biophys_palette = plot_data['PlotBiophys']['biophys_palette']

                if biophys_palette:
                    draw._biophys_palette = biophys_palette
                draw.plot_biophys(path=biophys_tpath)
            else:
                draw.plot_biophys()
            
            print('Physicochemical Composition of Interface Layers plot is done!\n')
        
        if 'PlotSASA' in plot_commands:

            if p_json:
                sasa_tpath = plot_data['PlotSASA']['sasa_path']
                draw.plot_SASA(path=sasa_tpath)

            else:
                draw.plot_SASA()
            
            print('Dynamic Evolution of Interface Area (SASA) plot is done!\n')

        if 'PlotiRMSD' in plot_commands:

            if p_json:
                irmsd_tpath = plot_data['PlotiRMSD']['irmsd_path']
                draw.plot_irmsd(path=irmsd_tpath)

            else:
                draw.plot_irmsd()

            print('Time Evolution of Interface Stability (i-RMSD) plot is done!\n')

        if 'PlotlRMSD' in plot_commands:

            if p_json:
                lrmsd_tpath = plot_data['PlotlRMSD']['lrmsd_path']
                draw.plot_lrmsd(path=lrmsd_tpath)

            else:
                draw.plot_lrmsd()

            print('Ligand RMSD (l-RMSD) Displacement Profile plot is done!\n')

        if 'PlotDockQ' in plot_commands:

            if p_json:
                dockq_tpath = plot_data['PlotDockQ']['dockq_path']
                draw.plot_dockq(path=dockq_tpath)

            else:
                draw.plot_dockq()

            print('Time Evolution of DockQ Score plot is done!\n')

        if 'PlotFnonnat' in plot_commands:

            if p_json:
                fnonnat_tpath = plot_data['PlotFnonnat']['fnonnat_path']
                draw.plot_fnonnat(path=fnonnat_tpath)

            else:
                draw.plot_fnonnat()

            print('Fraction of Non-Native Contacts (f_nonnat) plot is done!\n')
        
        if 'PlotFnat' in plot_commands:

            if p_json:
                fnat_tpath = plot_data['PlotFnat']['fnat_path']
                draw.plot_fnat(path=fnat_tpath)

            else:
                draw.plot_fnat()

            print('Fraction of Native Contacts (f_nat) plot is done!\n')

        if 'PlotPairwiseFreq' in plot_commands:

            if p_json:
                bondfreq_tpath = plot_data['PlotPairwiseFreq']['bar_table_path']
                bondfreq_palette = plot_data['PlotPairwiseFreq']['bar_palette']

                if bondfreq_palette:
                    draw._bar_palette = bondfreq_palette
                
                draw.plot_pairwise_freq(path=bondfreq_tpath)

            else:
                draw.plot_pairwise_freq()

            print('Frequency of Intermolecular Interaction plots are done!\n')

        if 'PlotResEne' in plot_commands:
            int_ene_thr = 50.0
            int_ene_itpath = None
            int_ene_tpath  = None

            if p_json:
                pe = plot_data.get('PlotResEne', {})
                int_ene_thr    = pe.get('interface_th', int_ene_thr)
                int_ene_itpath = pe.get('interface_table_path')
                int_ene_tpath  = pe.get('residue_based_table')
            
            success = draw.plot_int_energy(threshold=int_ene_thr, intf_path=int_ene_itpath, res_path=int_ene_tpath)
            if success:
                print('Distribution of Interface Residue Energies plots are done!\n')

        if 'PlotDSSP' in plot_commands:
            dssp_path = None
            dssp_threshold = 50.0
            dssp_intf_path = None

            if p_json:
                dssp_path = plot_data['PlotDSSP']['dssp_path']
                dssp_threshold = plot_data['PlotDSSP']['dssp_threshold']
                dssp_intf_path = plot_data['PlotDSSP']['dssp_intf_path']

            draw.plot_DSSP(path=dssp_path, threshold=dssp_threshold, intf_path=dssp_intf_path)
            print('Secondary Structure Stability at the Protein Interface plot is done!\n')


        if not p_json:
            draw._get_params_()
            print_stars(1)
            print("\n")
            print("Your plot run parameters are downloaded as plot_params.json")
            print("\n")
            print_stars(1)
            print("\n")
        else:
            f2.close()

if __name__ == '__main__':
    main()
