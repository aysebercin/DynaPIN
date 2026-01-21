import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from DynaPIN.handling import tables_errors
import json
from matplotlib import cm
from matplotlib.lines import Line2D
import warnings

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


class Plotter:
    def __init__(self, output_dir=None):
        """A class to visualize several analyses from dynapin class. Reads csv files and performs visualizations. Saves figures into the folder named 'figures', with 300 dpi and .png extension.
        
        Return: return_description
        """
        self.handler = tables_errors()

        #params
        self._output_dir = output_dir
        self._rmsd = False
        self._rmsd_path = None
        self._rg = False
        self._rg_path = None
        self._rmsf = False
        self._rmsf_path = None
        self._rmsf_intf_path = None
        self._biophys = False
        self._biophys_path = None
        self._biophys_palette = None
        self._pie = False
        self._pie_hbond = None
        self._pie_hph = None
        self._pie_ionic = None
        self._pie_path = None
        self._pie_palette = None
        self._bar = False
        self._bar_path = None
        self._bar_palette = None
        self._ene = False
        self._ene_thr = None
        self._ene_intf = None
        self._ene_path = None
        self._plot_SASA = False
        self._sasa_path = None
        self._plot_irmsd = False
        self._irmsd_path = None
        self._plot_fnonnat = False
        self._fnonnat_path = None
        self._plot_lrmsd = False
        self._lrmsd_path = None
        self._plot_dockq = False
        self._dockq_path = None
        self._plot_fnat = False
        self._fnat_path = None
        self._plot_dssp = False
        self._dssp_path = None
        self._dssp_intf_path = None
        self._dssp_threshold = 50.0

        sns.set_theme(context='notebook', style='whitegrid')
        plt.rcParams.update({
            'grid.linestyle': ':',
            'grid.linewidth': 0.6,
            'grid.color': '#d0d0d0',
            'axes.titleweight': 'bold',
            'axes.labelweight': 'bold',
            'figure.titleweight': 'bold',
            'axes.titlesize': 14,
            'axes.labelsize': 12,
            'xtick.labelsize': 10,
            'ytick.labelsize': 10,
            'legend.fontsize': 10,
        })

        if output_dir is None:
            output_dir = input('Plase provide job name, or existing job file where analysis folder is located:\n')
            self.job_path = output_dir

        else:
            self.job_path = os.path.abspath(output_dir)

        self.target_path = os.path.join(self.job_path, "figures")
        self.table_path = os.path.join(self.job_path, "tables")
        if not os.path.exists(self.target_path):
            os.mkdir(self.target_path)

        self.chain_colors = ['#648FFF', '#FFB000', '#DC267F', '#785EF0']
        self.complex_color = '#332288'

    def _get_chain_color(self, index):
        """Return a chain color, repeating the palette if needed."""
        if not self.chain_colors:
            return '#333333'
        return self.chain_colors[index % len(self.chain_colors)]

    def _category_palette(self, categories, palette_name="colorblind"):
        """Generate a color palette keyed by unique category names."""
        unique = list(dict.fromkeys(categories))
        if not unique:
            return {}
        colors = sns.color_palette(palette_name, len(unique))
        return {cat: colors[i % len(colors)] for i, cat in enumerate(unique)}

    def _place_axis_legend(self, ax, title=None, bbox=(1.02, 0.5), loc='center left'):
        """Place the legend outside the given axis."""
        handles, labels = ax.get_legend_handles_labels()
        if not handles:
            return None
        legend = ax.legend(handles, labels, loc=loc, bbox_to_anchor=bbox,
                           frameon=False, borderaxespad=0, title=title)
        return legend

    def _load_csv_safe(self, path, description):
        """Load CSV with graceful fallback; return None on error/empty."""
        try:
            self.handler.test_inp_path(path)
            df = pd.read_csv(path)
        except Exception as e:
            print(f"[warn] Skipping plot: could not read {description} ({e}).")
            return None
        if df.empty:
            print(f"[warn] Skipping plot: {description} is empty ({path}).")
            return None
        return df

    def plot_rmsd(self, path=None):
        """A function to perform line plot RMSD visualization for each chain and the overall complex. 'Reads qc_trajectory_metrics.csv' file.
        
        Return: None
        """
        
        
        if path:
            df = self._load_csv_safe(path, "RMSD table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "qc_trajectory_metrics.csv"), "RMSD table")
        if df is None:
            return

        self._rmsd = True
        self._rmsd_path = path

        time_name = df.columns[0]
        time_var = df.iloc[:, 0].values


        cols_use = [x for x in df.columns.tolist() if "RMSD" in x and len(x) > 5]
        chains = cols_use[:-1]
        complex = cols_use[-1]

        fig, ax = plt.subplots(figsize=(8, 4.5), layout='constrained')
        for i, col in enumerate(chains):
            mol_name = col.split("RMSD")
            ax.plot(time_var, df[col], label=mol_name[0], color=self._get_chain_color(i))
        
        ax.plot(time_var, df[complex], label='Complex', color=self.complex_color)
        ax.set_xlabel(time_name)
        ax.set_ylabel(f'RMSD (Å)')

        ax.set_title(f"Time Evolution of backbone Root Mean Square Deviation (RMSD)")
        self._place_axis_legend(ax)
        fig.savefig(os.path.join(self.target_path, f'qc_rmsd_backbone.png'), dpi=300, bbox_inches='tight')

    def plot_rg(self, path=None):
        """A function to perform line plot RG visualization for each chain and the overall complex. 'Reads qc_trajectory_metrics.csv' file.
        
        Return: None
        """
        if path:
            df = self._load_csv_safe(path, "RG table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "qc_trajectory_metrics.csv"), "RG table")
        if df is None:
            return

        self._rg = True
        self._rg_path = path

        cols_use = [x for x in df.columns.tolist() if "RG" in x]
        chains = cols_use[:-1]
        complex = cols_use[-1]

        time_name = df.columns[0]
        time_var = df.iloc[:, 0].values

        fig, ax = plt.subplots(figsize=(8, 4.5), layout='constrained')
        for i, col in enumerate(chains):
            mol_name = col.split("RG")
            ax.plot(time_var, df[col], label=mol_name[0], color=self._get_chain_color(i))

        ax.plot(time_var, df[complex], label='Complex', color=self.complex_color)
        ax.set_xlabel(time_name)
        ax.set_ylabel(f'$R_g$ (Å)')

        ax.set_title(f"Time Evolution of Radius of Gyration ($R_g$)")
        self._place_axis_legend(ax)
        fig.savefig(os.path.join(self.target_path, f'qc_rg_complex.png'), dpi=300, bbox_inches='tight')

    def plot_rmsf(self, rmsf_path=None, intf_path=None):
        """
        A function to perform line plot RMSF visualization for each chain, with interface residues marked.
        Reads 'qc_residue_rmsf.csv' and optionally 'res_interface_stats.csv'.
        """
        rmsf_df_path = rmsf_path or os.path.join(self.table_path, "qc_residue_rmsf.csv")
        df1 = self._load_csv_safe(rmsf_df_path, "RMSF table")
        if df1 is None:
            return

        g2 = None
        try:
            path_to_load = None
            if intf_path:
                if os.path.exists(intf_path):
                    path_to_load = intf_path
                else:
                    warnings.warn(f"Provided interface file not found: {intf_path}. Plotting without markers.")
            else:
                default_path = os.path.join(self.table_path, "res_interface_stats.csv")
                if os.path.exists(default_path):
                    path_to_load = default_path
            
            if path_to_load:
                intf_df = self._load_interface_table(path_to_load)
                if intf_df is not None and not intf_df.empty:
                    mask = (intf_df["Interface Label"].isin([2, 3, 4])) & (intf_df["Interface_score"] >= 50.0)
                    filtered_intf_df = intf_df[mask]
                    if not filtered_intf_df.empty:
                        g2 = filtered_intf_df.groupby("Chain")
                    else:
                        warnings.warn("No residues in the interface file met the filter criteria. Plotting without markers.")
                else:
                     warnings.warn("Interface file is empty or invalid. Plotting without markers.")
        except Exception as e:
            warnings.warn(f"An error occurred while processing interface data: {e}. Plotting without markers.")

        self._rmsf = True
        self._rmsf_path = rmsf_path
        self._rmsf_intf_path = intf_path

        groups = df1.groupby("Molecule")
        molecule_names = np.unique(df1["Molecule"])
        n_groups = len(molecule_names)
        residue_counts = [len(groups.get_group(name)) for name in molecule_names]
        
        fig_width = np.clip(sum(residue_counts) * 0.04, 8, 30)
        fig, axes = plt.subplots(1, n_groups, figsize=(fig_width, 6), sharey=True,
                                 gridspec_kw={'width_ratios': residue_counts})
        if n_groups == 1:
            axes = [axes]

        any_interface_markers = False
        for ind, (ax, chain_name) in enumerate(zip(axes, molecule_names)):
            chain_data = groups.get_group(chain_name)
            chain_color = self._get_chain_color(ind)
            
            interface_res_nums = set()
            if g2 is not None:
                try:
                    int_chain_data = g2.get_group(chain_name)
                    interface_res_nums = {int(res[3:]) for res in int_chain_data["Residue"]}
                except KeyError:
                    pass  

            # Plot the main RMSF line
            ax.plot(chain_data["Residue Number"], chain_data["RMSF"], color=chain_color, linewidth=2)

            # Find which data points are interface residues and plot markers
            if interface_res_nums:
                marker_indices = chain_data["Residue Number"].isin(interface_res_nums)
                interface_points = chain_data[marker_indices]
                if not interface_points.empty:
                    ax.plot(interface_points["Residue Number"], interface_points["RMSF"], 'o', 
                            markerfacecolor="red", markeredgecolor="red", 
                            linestyle='None', markersize=4)
                    any_interface_markers = True

            ax.set_xlabel("Residue Number")
            ax.set_title(f"Chain {chain_name}")

        axes[0].set_ylabel("RMSF (Å)")
        plt.suptitle("Residue-based backbone Root Mean Square Fluctuation (RMSF)")
        
        if any_interface_markers:
            legend_ax = axes[1] if n_groups > 1 else axes[0]
            interface_handle = Line2D([], [], color='red', marker='o', linestyle='None',
                                      markersize=5, label='Dynamic Interface Residues')
            legend_ax.legend(handles=[interface_handle], loc='upper right')

        fig.savefig(os.path.join(self.target_path, 'qc_rmsf_backbone.png'), dpi=300, bbox_inches='tight')
        plt.close(fig)

    def plot_irmsd(self, path=None):
        """A function to perform lineplot visualization of interface RMSD through the simulation. Reads 'qc_trajectory_metrics.csv' file.
        
        Return: None
        """

        self._plot_irmsd = True

        if path:
            df = self._load_csv_safe(path, "iRMSD table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "qc_trajectory_metrics.csv"), "iRMSD table")
        if df is None:
            return

        
        self._irmsd_path = path
        
        if "DockQ_mapping" not in df.columns:
            print("[warn] Skipping iRMSD plot; 'DockQ_mapping' column not found.")
            print("      Hint: This is expected for monomer systems.")
            return

        d = df.groupby("DockQ_mapping")
        palette = self._category_palette(d.groups.keys())

        fig,ax  = plt.subplots(figsize=(8, 4.5), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            color = palette.get(i, self.complex_color)
            ax.plot(data.columns[0], "iRMSD", data=data, label=i, color=self.complex_color)
        
        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('$i$-RMSD (Å)')
        ax.set_ylim(bottom=0)

        ax.set_title('Time Evolution of Interface Stability ($i$-RMSD)')
        self._place_axis_legend(ax, title="Mapping")
        fig.savefig(os.path.join(self.target_path, f'qc_irmsd_interface.png'), dpi=300, bbox_inches='tight')

    def plot_lrmsd(self, path=None):
        """A function to perform lineplot visualization of ligand RMSD through the simulation. Reads 'qc_trajectory_metrics.csv' file.
        
        Return: None
        """

        self._plot_lrmsd = True

        if path:
            df = self._load_csv_safe(path, "lRMSD table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "qc_trajectory_metrics.csv"), "lRMSD table")
        if df is None:
            return

        
        self._lrmsd_path = path

        if "DockQ_mapping" not in df.columns:
            print("[warn] Skipping lRMSD plot; 'DockQ_mapping' column not found.")
            print("      Hint: This is expected for monomer systems.")
            return

        d = df.groupby("DockQ_mapping")
        palette = self._category_palette(d.groups.keys())

        fig,ax  = plt.subplots(figsize=(8, 4.5), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            color = palette.get(i, self.complex_color)
            ax.plot(data.columns[0], "lRMSD", data=data, label=i, color=self.complex_color)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('$l$-RMSD (Å)')
        ax.set_ylim(bottom=0)


        ax.set_title('Time Evolution of Ligand Displacement Profile ($l$-RMSD)')
        self._place_axis_legend(ax, title="Mapping")
        fig.savefig(os.path.join(self.target_path, f'qc_lrmsd_ligand.png'), dpi=300, bbox_inches='tight')

    def plot_dockq(self, path=None):
        """A function to perform lineplot visualization of ligand RMSD through the simulation. Reads 'qc_trajectory_metrics.csv' file.
        
        Return: None
        """

        self._plot_dockq = True

        if path:
            df = self._load_csv_safe(path, "DockQ table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "qc_trajectory_metrics.csv"), "DockQ table")
        if df is None:
            return

        
        self._dockq_path = path

        if "DockQ_mapping" not in df.columns:
            print("[warn] Skipping DockQ plot; 'DockQ_mapping' column not found.")
            print("      Hint: This is expected for monomer systems.")
            return

        d = df.groupby("DockQ_mapping")
        palette = self._category_palette(d.groups.keys())

        fig,ax  = plt.subplots(figsize=(8, 4.5), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            color = palette.get(i, self.complex_color)
            ax.plot(data.columns[0], "DockQ", data=data, label=i, color=self.complex_color)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('DockQ')
        ax.set_ylim(bottom=0)

        self._place_axis_legend(ax, title="Mapping")
        ax.set_title('Time Evolution of DockQ Score')
        fig.savefig(os.path.join(self.target_path, f'qc_dockq_score.png'), dpi=300, bbox_inches='tight')

    def plot_fnonnat(self, path=None):
        """A function to perform lineplot visualization of fraction of native contacts through the simulation. Reads 'qc_trajectory_metrics.csv' file.
        
        Return: None
        """

        self._plot_fnonnat = True

        if path:
            df = self._load_csv_safe(path, "Fnonnat table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "qc_trajectory_metrics.csv"), "Fnonnat table")
        if df is None:
            return

        self._fnonnat_path = path
        
        if "DockQ_mapping" not in df.columns:
            print("[warn] Skipping FNONNAT plot; 'DockQ_mapping' column not found.")
            print("      Hint: This is expected for monomer systems.")
            return

        d = df.groupby("DockQ_mapping")
        palette = self._category_palette(d.groups.keys())

        fig,ax  = plt.subplots(figsize=(8, 4.5), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            color = palette.get(i, self.complex_color)
            ax.plot(data.columns[0], "fnonnat", data=data, label=i, color=self.complex_color)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('$f_{nonnat}$')
        ax.set_ylim(bottom=0)

        self._place_axis_legend(ax, title="Mapping")
        ax.set_title('Fraction of Non-Native Contacts ($f_{nonnat}$)')
        fig.savefig(os.path.join(self.target_path, f'qc_fnonnat_contacts.png'), dpi=300, bbox_inches='tight')

    def plot_fnat(self, path=None):
        """A function to perform lineplot visualization of fraction of native contacts through the simulation. Reads 'qc_trajectory_metrics.csv' file.
        
        Return: None
        """

        self._plot_fnat = True

        if path:
            df = self._load_csv_safe(path, "Fnat table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "qc_trajectory_metrics.csv"), "Fnat table")
        if df is None:
            return

        self._fnat_path = path

        if "DockQ_mapping" not in df.columns:
            print("[warn] Skipping FNAT plot; 'DockQ_mapping' column not found.")
            print("      Hint: This is expected for monomer systems.")
            return

        d = df.groupby("DockQ_mapping")
        palette = self._category_palette(d.groups.keys())

        fig,ax  = plt.subplots(figsize=(8, 4.5), layout='constrained')

        for i in d.groups:
            data = d.get_group(i)
            color = palette.get(i, self.complex_color)
            ax.plot(data.columns[0], "fnat", data=data, label=i, color=self.complex_color)

        ax.set_xlabel(df.columns[0])
        ax.set_ylabel('$f_{nat}$')
        ax.set_ylim(bottom=0)

        self._place_axis_legend(ax, title="Mapping")
        ax.set_title('Fraction of Native Contacts ($f_{nat}$)')
        fig.savefig(os.path.join(self.target_path, f'qc_fnat_contacts.png'), dpi=300, bbox_inches='tight')

    def plot_biophys(self, path=None):
        """A function to perform barplot core-rim and biophysical classification visualization. Reads 'inetrface_label.csv' file.
        
        Return: None
        """
        if not self._biophys_palette:
            self._biophys_palette = sns.color_palette("Dark2", 6)
        

        if path:
            intf_df = self._load_interface_table(path)
        else:    
            intf_df = self._load_interface_table()
        if intf_df is None or intf_df.empty:
            print("[warn] Skipping Biophysical plot: interface table missing or empty.")
            return

        self._biophys = True
        self._biophys_path = path

        def _calc_biophys_type(res_type):
            res_dict = {"GLY": 'hydrophobic', "ALA": 'hydrophobic', "PRO": 'hydrophobic', "VAL": 'hydrophobic',
                        "LEU": 'hydrophobic', "ILE": 'hydrophobic', "MET": 'hydrophobic', "TRP": 'hydrophobic',
                        "PHE": 'hydrophobic',
                        "SER": 'polar',
                        "THR": 'polar', "TYR": 'polar', "ASN": 'polar', "GLN": 'polar', "CYS": 'polar',
                        "LYS": '+ly charged', "ARG": '+ly charged', "HIS": '+ly charged', "ASP": '-ly charged',
                        "GLU": '-ly charged'}
            return res_dict.get(res_type, 'Other')

        biophy = list()
        for index, row in intf_df.iterrows():
            biophy.append(_calc_biophys_type(row["Residue"][:3]))

        intf_df["Biophysical Type"] = biophy

        plot_df = (
            intf_df.groupby(["Biophysical Type", "Interface Label"])
                   .size()
                   .reset_index(name="Count")
        )
        ty = list()
        for ind, r in plot_df.iterrows():
            if r["Interface Label"] == 2:
                ty.append("Support")
            elif r["Interface Label"] == 3:
                ty.append("Rim")
            elif r["Interface Label"] == 4:
                ty.append("Core")
            elif r["Interface Label"] == 1:
                ty.append("Surface")
            elif r["Interface Label"] == 0:
                ty.append("Interior")

        plot_df["Interface Text"] = ty

        total = plot_df["Count"].sum()
        plot_df["Percentage"] = plot_df["Count"] / total * 100
        plot_df.sort_values("Interface Label", inplace=True)

        num_bars = plot_df["Interface Text"].nunique()
        fig_width = max(8.0, num_bars * 1.0)
        fig, ax = plt.subplots(figsize=(fig_width, 6))
        sns.barplot(x="Interface Text", data=plot_df, hue="Biophysical Type", y="Percentage", ax=ax, palette=sns.color_palette(self._biophys_palette, 4))
        l = ["Support", "Rim", "Core"]
        for ind, i in enumerate(ax.get_xticklabels()):
            if ax.get_xticklabels()[ind].get_text() in l:
                plt.setp(ax.get_xticklabels()[ind], weight='bold')

        plt.ylabel("Occurance Rate of Residue Types (%)")
        plt.xlabel("Interface Label")

        plt.title("Physicochemical Composition of Interface Layers (Levy, 2010)")
        self._place_axis_legend(ax, title="Biophysical Type", bbox=(1.02, 1), loc='upper left')
        fig.savefig(os.path.join(self.target_path, f'res_biophys_composition.png'), dpi=300, bbox_inches='tight')

    def plot_DSSP(self, threshold=50.0, intf_path=None, path=None):
        self._plot_DSSP = True 
        self._dssp_threshold=threshold
        
        if path:
            df = self._load_csv_safe(path, "DSSP residue table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "res_trajectory_props.csv"), "DSSP residue table")
        if df is None:
            return

        if intf_path:
            intf_df = self._load_interface_table(intf_path)
        else:
            intf_df = self._load_interface_table()
        if intf_df is None or intf_df.empty:
            print("[warn] Interface table missing/empty for DSSP plot; skipping.")
            return

        self._dssp_path = path  
        self._dssp_intf_path = intf_path

        time_name = df.columns[0]

        mask = (intf_df["Interface Label"].isin([2, 3, 4])) & (intf_df["Interface_score"] >= float(threshold))
        int_df = intf_df[mask]

        df = df.loc[:, [time_name,"Chain",'Residue', 'Residue Number', 'Secondary Structure']]

        df["Residue"] = [a + str(b) for a, b in zip(df["Residue"], df["Residue Number"])]
        mylabels = {
            'H': 'Alpha Helix',
            'B': 'Beta Bridge',
            'E': 'Strand',
            'G': 'Helix-3',
            'I': 'Helix-5',
            'T': 'Turn',
            'S': 'Bend',
            'P': 'Coil',
            'C': 'Coil',
            'Null': 'Null',
            '!': 'Chain Break'
        }

        dssp_palette = {
            'H': '#648FFF',
            'B': '#FFB000',
            'E': '#DC267F',
            'G': '#785EF0',
            'I': '#44AA99',
            'T': '#FE6100',
            'S': '#e377c2',
            'P': '#7f7f7f',
            'C': '#bcbd22',
            'Null': '#17becf',
            '!': '#000000'
        }

        # Load data
        interface_groupped = int_df.groupby('Chain')

        groups = df.groupby(['Chain'])
        n_groups = len(groups)
        fig, axes = plt.subplots(n_groups, 1, figsize=(15, 15))
        if n_groups == 1:
            axes = [axes]
        used_codes = set()
        
        for ax, chain, interface in zip(axes, groups, interface_groupped):
        
            chain_df = chain[1]

            intf_ch = interface[1]

            intf_residues = intf_ch['Residue'].tolist()
            chain_df = chain_df[chain_df['Residue'].isin(intf_residues)]

            sns.scatterplot(
                data=chain_df,
                x=chain_df.columns[0],
                y="Residue",
                hue="Secondary Structure",
                palette=dssp_palette,
                marker="+",
                linewidths=2,
                ax=ax,
                s=50
            )
            used_codes.update(chain_df["Secondary Structure"].unique().tolist())

            ax.set_xlabel(f'{time_name}')
            ax.set_ylabel('Residue')
            ax.set_title(f"Chain {chain[0][0]}")
            ax.get_legend().remove()

        plt.suptitle("Secondary Structure Profile of Interface Residues")
        legend_handles = []
        legend_labels = []
        for code, label in mylabels.items():
            if code in used_codes:
                color = dssp_palette.get(code, '#888888')
                legend_handles.append(Line2D([0], [0], marker='s', linestyle='', markersize=8,
                                              markerfacecolor=color, markeredgecolor=color))
                legend_labels.append(label)
        if legend_handles:
            fig.legend(legend_handles, legend_labels, ncol=1, loc='upper left', bbox_to_anchor=(1.02, 1),
                       frameon=False, title="Secondary Structure")
        plt.tight_layout()
        fig.subplots_adjust(hspace=0.20)
        fig.savefig(os.path.join(self.target_path, f'res_dssp_stability.png'), dpi=300, bbox_inches='tight')

    def plot_SASA(self, path=None):
        """A function to perform line plot Interface Area in A^2 visualization for each chain and the overall complex. 'Reads res_trajectory_props.csv' file.
        
        Return: None
        """

        self.plot_SASA = True

        if path:
            df = self._load_csv_safe(path, "SASA table")
        else:
            df = self._load_csv_safe(os.path.join(self.table_path, "res_trajectory_props.csv"), "SASA table")
        if df is None:
            return

        self._sasa_path = path
        df = df.iloc[:, [0,1,2,3,7,8]]

        time_name = df.columns[0]

        plot_df = pd.DataFrame(columns=[time_name, "Chain", "SASA"])


        int_df = df[(df["Interface Label"] == 4) | (df["Interface Label"] == 2) | (df["Interface Label"] == 3)]
            
        g2 = int_df.groupby(["Chain"])

        fig, ax = plt.subplots(figsize=(8, 4.5), layout='constrained')
        for chain in g2.groups:
            chain_df = g2.get_group(chain)
            frames = chain_df.groupby([time_name])
            
            for frame in frames.groups:
                frame_df = frames.get_group(frame)
                frame_sasa = frame_df['SASA'].sum()
                plot_df.loc[len(plot_df)] = [frame, chain, frame_sasa]

        for i, chain_ in enumerate(plot_df.groupby('Chain').groups):
            ax.plot(plot_df.groupby('Chain').get_group(chain_)[time_name], plot_df.groupby('Chain').get_group(chain_)['SASA'], label=chain_, color=self._get_chain_color(i))

        sums = [plot_df.groupby(time_name).get_group(group)['SASA'].sum() for group in plot_df.groupby(time_name).groups]

        ax.plot(plot_df[time_name].unique() ,sums, color=self.complex_color, label="Complex")
         
        ax.set_xlabel(time_name)
        ax.set_ylabel(f'Interface Area (Å²)')

        ax.set_title(f"Time Evolution of Solvent Accessible Surface Area of Interface")
        self._place_axis_legend(ax, title="Chains")
        fig.savefig(os.path.join(self.target_path, f'res_sasa_interface.png'), dpi=300, bbox_inches='tight')

    def plot_pairwise_freq(self, path=None):
        """A function to perform barplot visualization of interaction frequency in simulation of residue pairs with locations (side chain or backbone) of atoms that participated in interaction. Reads 'int_pairwise_trajectory.csv' file.
        
        Return: None
        """
        interaction_titles = {
            'hbond': 'hydrogen bonds',
            'hydrophobic': 'hydrophobic interactions',
            'ionic': 'salt bridges'
        }

        if path:
            base_df = self._load_csv_safe(path, "Interaction table")
        else:
            base_df = self._load_csv_safe(os.path.join(self.table_path, "int_pairwise_trajectory.csv"), "Interaction table")
        if base_df is None:
            return
        if base_df.empty:
            print("[warn] Interaction table is empty; skipping pairwise frequency plot.")
            return
        time_name = base_df.columns[0]
        time_var = base_df.iloc[:,0]
        sim_time = len(np.unique(time_var.tolist()))
        df = base_df[[base_df.columns[0],"itype", "pairwise", "atom_a", "atom_b"]].groupby('itype')
            
        self._bar = True
        self._bar_path = path

        if not self._bar_palette:
            self._bar_palette = [
                            # 1) H-bonds
                            ["#A2C9EF",  # light  – bb–bb
                            "#5286BA",  # medium – bb–sc / sc–bb
                            "#1C3D5E"], # darker – sc–sc

                            # 2) Hydrophobic
                            ["#ACF4E8",  # light  – bb–bb
                            "#70C3B4",  # medium – bb–sc / sc–bb
                            "#348D8C"], # darker – sc–sc

                            # 3) Ionic
                            ["#F3ABCF",  # light  – bb–bb
                            "#CC6677",  # medium – bb–sc / sc–bb
                            "#882255"], # darker – sc–sc
                        ]
                             
        backbones = {"HN", "N", "CA", "HA", "C", "O"}

        for x,palet in zip(df.groups, self._bar_palette):
            ll = dict()
            z = df.get_group(x).copy()
            z_gr = z.groupby('pairwise')
            for pair in z_gr.groups:
                df2 = z_gr.get_group(pair)
                l = len(np.unique(df2.iloc[:,0].tolist()))
                percentage = l/sim_time*100
                ll[pair] = percentage

            def is_backbone(atom_name):
                if atom_name is None:
                    return False
                cleaned = str(atom_name).strip().upper()
                return cleaned in backbones

            def classify_loc(atom_a, atom_b):
                a_bb = is_backbone(atom_a)
                b_bb = is_backbone(atom_b)
                if a_bb and b_bb:
                    return 0  # bb-bb
                if a_bb or b_bb:
                    return 1  # bb-sc or sc-bb
                return 2      # sc-sc

            plot_dict = dict()
            for _, row in z.iterrows():
                pair_key = row["pairwise"]
                if pair_key not in plot_dict:
                    plot_dict[pair_key] = [0, 0, 0]
                loc_idx = classify_loc(row["atom_a"], row["atom_b"])
                plot_dict[pair_key][loc_idx] += 1

            entries = []
            for pair, counts in plot_dict.items():
                total = sum(counts)
                if total == 0:
                    scaled = [0.0, 0.0, 0.0]
                else:
                    scaled = [(counts[i] * ll[pair]) / total for i in range(3)]
                entries.append((pair, scaled))

            entries.sort(key=lambda item: sum(item[1]), reverse=True)

            fixed_height = 8.0
            loc_labels = ("bb-bb", "bb-sc or sc-bb", "sc-sc")
            bar_width = 0.7
            interaction_key = x.lower().strip()
            full_name = interaction_titles.get(interaction_key, x.strip().capitalize())

            if not entries:
                continue

            chunk_pairs = [item[0] for item in entries]
            chunk_vals = np.array([item[1] for item in entries])

            fig_width = max(12.0, len(chunk_pairs) * 0.8)
            fig, ax = plt.subplots(figsize=(fig_width, fixed_height))

            bottom = np.zeros(len(chunk_pairs))
            colors = palet
            for idx, label in enumerate(loc_labels):
                ax.bar(chunk_pairs, chunk_vals[:, idx], width=bar_width,
                       bottom=bottom, color=colors[idx], label=label)
                bottom += chunk_vals[:, idx]

            self._place_axis_legend(ax, bbox=(1.01, 1), loc='upper left', title="Interaction Type")

            plt.xlabel("Interacting Pairs")
            plt.ylabel(f"Persistence rate of {full_name} (%)")
            plt.xticks(rotation=85)
            plt.yticks([tick for tick in range(0, 101, 10)])
            plt.ylim(ymin=0, ymax=105)

            plt.title(f"Frequency of Intermolecular {full_name} across simulation")

            outfile = os.path.join(self.target_path, f"int_pairwise_{x}_frequency.png")
            plt.savefig(outfile, dpi=300, bbox_inches='tight', format="png")
            plt.close(fig)

    def plot_int_energy(self, threshold=50.0, intf_path=None, res_path=None):
        """A function to perform boxplot visualization of interaction energy variation of residues that has constant interface label at least 50%-by default- of all simulation. Reads 'res_interface_stats' and 'res_trajectory_props.csv' files.
        
        Keyword arguments:
        threshold -- int: The desired percentage value of the residual interface label to be constant throughout the simulation
        Return: None
        """
        if not res_path:
            res_path = os.path.join(self.table_path, "res_trajectory_props.csv")

        intf_df = self._load_interface_table(intf_path)
        if intf_df is None or intf_df.empty:
            print("[warn] Skipping PlotResEne; interface table missing or empty.")
            return False
        energy_df_full = self._load_csv_safe(res_path, "residue-energy table")
        if energy_df_full is None:
            return False

        energy_cols = [c for c in ["Total Residue Energy", "Van der Waals Energy", "Electrostatic Energy"] if c in energy_df_full.columns]
        if not energy_cols:
            print("[warn] Skipping PlotResEne: no energy columns found (expected any of Total Residue Energy, Van der Waals Energy, Electrostatic Energy).")
            print("      Hint: Run the analysis with --foldx_path to generate residue energies of your complex.")
            return False

        energy_df = energy_df_full[["Chain", "Residue", "Residue Number"] + energy_cols].copy()

        self._ene = True
        self._ene_thr = threshold
        self._ene_intf = intf_path
        self._ene_path = res_path

        mask = (intf_df["Interface Label"].isin([2, 3, 4])) & (intf_df["Interface_score"] >= float(threshold))
        int_df_filtered = intf_df[mask]

        energy_df["Residue"] = [a + str(b) for a, b in zip(energy_df["Residue"], energy_df["Residue Number"])]
        g = energy_df.groupby(["Chain", "Residue"])
        n_df = pd.DataFrame(columns=["Chain", "Residue", "Residue Number"] + energy_cols)
        my_palette = dict()
        ch_num = 1
        for index, row in int_df_filtered.iterrows():
            mg = (row["Chain"], row["Residue"])
            if mg[0] not in my_palette.keys():
                my_palette[mg[0]] = self._get_chain_color(ch_num - 1)
                ch_num += 1
            try:
                data = g.get_group(mg)
                n_df = pd.concat([n_df, data])
            except KeyError:
                continue # Skip if residue from interface file is not in energy file

        if n_df.empty:
            print("[warn] No matching interface residues found in the energy data. Skipping PlotResEne.")
            return False
            
        plot_data = n_df.sort_values(["Chain", "Residue Number"])

        for col in energy_cols:
            fig, ax = plt.subplots(figsize=(15, 8))
            sns.boxplot(data=plot_data, x="Residue", y=col, hue="Chain", palette=my_palette)
            plt.xticks(rotation=90)
            energy_units = {
                "Van der Waals Energy": "kcal/mol",
                "Electrostatic Energy": "kcal/mol",
                "Total Residue Energy": "a.u",
            }
            if col in energy_units:
                ax.set_ylabel(f"{col} ({energy_units[col]})")
            plt.title(f"Distribution of Interface Residue Energies: {col} ")
            self._place_axis_legend(ax, title="Chain", bbox=(1.02, 1), loc='upper left')
            fname = col.lower().replace(" ", "_")
            plt.savefig(os.path.join(self.target_path, f"res_energy_{fname}.png"), dpi=300, bbox_inches='tight')
            plt.close(fig)
        return True

    def _load_interface_table(self, path=None):
        """
        Load res_interface_stats.csv and derive compatibility columns used by older plots.
        """
        csv_path = path if path else os.path.join(self.table_path, "res_interface_stats.csv")
        description = os.path.basename(csv_path)
        df = self._load_csv_safe(csv_path, description)
        if df is None:
            return None
        label_cols = [f"label_{i}" for i in range(5) if f"label_{i}" in df.columns]
        if not label_cols:
            print("[warn] res_interface_stats.csv is missing label columns (label_0 ... label_4). Skipping plot.")
            return None

        for col in label_cols:
            df[col] = pd.to_numeric(df[col], errors="coerce").fillna(0.0)

        if "interface" in df.columns:
            df["Interface_score"] = pd.to_numeric(df["interface"], errors="coerce").fillna(0.0)
        else:
            df["Interface_score"] = df[[c for c in label_cols if c in ("label_2", "label_3", "label_4")]].sum(axis=1)

        df["Interface Label"] = (
            df[label_cols]
              .astype(float)
              .idxmax(axis=1)
              .str.replace("label_", "", regex=False)
              .astype(int)
        )

        return df

    def _get_params_(self):
        params = {
            'output_dir':self._output_dir,
            'output_path': self.target_path,
            'table_path': self.table_path,
            'chain_color_palette': self.chain_colors,
            'complex_color': self.complex_color,
            'PlotRMSD': {
                'Run': self._rmsd,
                'rmsd_table_path':self._rmsd_path
                },
            'PlotRG': {
                'Run': self._rg,
                'rg_table_path':self._rg_path
                },
            'PlotRMSF': {
                'Run': self._rmsf,
                'rmsf_table_path':self._rmsf_path,
                'rmsf_intf_res_table': self._rmsf_intf_path
                },
            'PlotBiophys': {
                'Run': self._biophys,
                'biophys_table_path':self._biophys_path,
                'biophys_palette': self._biophys_palette
                },
            'PlotSASA' : {
                'Run': self._plot_SASA,
                'sasa_path': self._sasa_path,
            },
            'PlotiRMSD': {
                'Run': self._plot_irmsd,
                'irmsd_path': self._irmsd_path,

            },
            'PlotlRMSD': {
                'Run': self._plot_lrmsd,
                'lrmsd_path': self._lrmsd_path,

            },
            'PlotDockQ': {
                'Run': self._plot_dockq,
                'dockq_path': self._dockq_path,

            },
            'PlotFnonnat': {
                'Run': self._plot_fnonnat,
                'fnonnat_path': self._fnonnat_path
            },
            'PlotFnat': {
                'Run': self._plot_fnat,
                'fnat_path': self._fnat_path
            },
            'PlotPairwiseFreq': {
                'Run': self._bar,
                'bar_table_path':self._bar_path,
                'bar_palette': self._bar_palette
                },
            'PlotResEne': {
                'Run': self._ene,
                'interface_th':self._ene_thr,
                'interface_table_path': self._ene_intf,
                'residue_based_table': self._ene_path
                },
            'PlotDSSP': {
                'Run': self._plot_dssp,
                'dssp_path': self._dssp_path,
                'dssp_threshold':self._dssp_threshold,
                'dssp_intf_path':self._dssp_intf_path
            },

            }

        for group_name, members in PLOT_GROUPS.items():
            params[group_name] = {
                'Run': any(params.get(member, {}).get('Run') for member in members),
                'members': members
            }

        json_path = os.path.join(self.job_path, 'plot_params.json')
        with open(json_path, 'w+') as ofh:
            json.dump(params, ofh)
