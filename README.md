# MetworkPy Application Note

Data and Code associated with 'MetworkPy: A Python Package for Graph- and
Information-theoretic Investigation of Metabolic Networks'.

## Usage

In order to run the analysis of the simulation model and the *Mycobacterium
tuberculosis*(Mtb) transcription factors (TF), first install the required
dependencies, see the [Dependencies section](#dependencies) below.

If you have pixi installed, the analysis pipelines can be run using pixi tasks.
Specifically,

```{bash}
# Run the analysis for the simulation model
pixi run simulation # Relatively Short
# Run the analysis for the Mtb TFs
pixi run mtb_tf # Very Long
```

Note that the simulation model analysis should take less than 10 minutes, but
the Mtb TF analysis will take much longer (on the order of 24 hours).

The analysis can be configured by modifying the values in the CONFIG.toml file
(or by directly modifying the python scripts in the scripts directory).

For alternative methods of installing dependencies, make sure you activate the
environment (either conda or venv), and then the scripts can be run directly by
calling

```{bash}
python <script>
```

in the script directory (there are some relative imports, so python has to be
run from with directory that directly contains the script of interest).

The scripts are numbered, and should be executed in order since some later ones
depend on earlier ones.

## Dependencies

This project mainly uses [pixi](https://pixi.prefix.dev/latest/) to manage
dependencies, but includes several different alternatives for installing
dependencies, including [conda](#conda-environment), [pip](#pip), and
[pixi](#pixi).

### Conda Environment

Conda can be installed using [Anaconda](https://www.anaconda.com/download),
[Miniconda](https://www.anaconda.com/docs/getting-started/miniconda/main), or
[Miniforge](https://github.com/conda-forge/miniforge). We recommend using
Miniforge as that is what this repository was tested with, but any method should
work.

After installing conda using one of the above methods, then the environment can
be created using the environment.yml file:

```{bash}
conda create --file environment.yml
```

and then activated with

```{bash}
conda activate metworkpy_application_note
```

### Pip

Make sure python is installed, you could use:

- [python download page](https://www.python.org/downloads/)
- [pyenv](https://github.com/pyenv/pyenv)
- [mise-en-place](https://mise.jdx.dev/lang/python.html)
- [brew](https://brew.sh/)
- [scoop](https://scoop.sh/)

or other package manager to install python. Alternatively,
[uv](https://docs.astral.sh/uv/) could be used to install a python version or to
directly create a virtual environment.

Once python is installed, a virtual environment can be created by

```{bash}
python -m venv .venv
```

then activated with

```{bash}
# Windows
.venv\Scripts\activate
# Macos/Linux
source myfirstproject/bin/activate
```

and finally required dependencies can be installed with

```{bash}
pip install -r requirements.txt
```

### Pixi

Pixi can be installed based on instructions from
[https://pixi.prefix.dev/latest/installation/](https://pixi.prefix.dev/latest/installation/).

## Directory Structure

<details>
<summary>Descriptions of the files in this repository</summary>
- data: Directory containing data for performing the *Mycobacterium
  tuberculosis* (Mtb) analysis, see [below](#data-references) for the full
  references.
  - bigg_info: Information from the [BiGG database](http://bigg.ucsd.edu/)
    - bigg_models_metabolites.txt: Information about metabolites in the BiGG
      database, used for contextualizing predictions about metabolites, and
      translating between iEK1011_v2 and iEK1011 reaction/metabolite ids.
  - gene_info: Information about Mtb genes
    - Mycobacterium_tuberculosis_H37Rv_txt_v5.txt: Information about
      mycobacterium tuberculosis genes, found at
      [https://mycobrowser.epfl.ch/](https://mycobrowser.epfl.ch/), used for
      gene information
    - bosch_vi.xlsx: Information about *Mycobacterium tuberculosis* gene
      essentiality and vulnerability index from
      [Genome-wide gene expression tuning reveals diverse vulnerabilities of M. tuberculosis](https://pmc.ncbi.nlm.nih.gov/articles/PMC8382161/)
  - mtb_transcription_factors: Information about the Transcription Factors (TF)
    in Mtb
    - tfoe_targets.xlsx: Information on the gene regulatory targets of the TFs
      in Mtb, from
      [Mapping and manipulating the Mycobacterium tuberculosis transcriptome using a transcription factor overexpression-derived regulatory network](https://pmc.ncbi.nlm.nih.gov/articles/PMC4249609/)
      used for determining the Mtb TF gene regulatory targets
- escher_maps: Base [Escher maps](https://escher.github.io/) for the simulation
  model and the iEK1011 Mtb model
- models: Genome Scale Metabolic Models (GSMM) for Mtb and an example simulation
  model. The Mtb models (iEK1011_v2) are from
  [A systematic evaluation of Mycobacterium tuberculosis Genome-Scale Metabolic Networks](https://pmc.ncbi.nlm.nih.gov/articles/PMC7316355/),
  modified for different media conditions (7H9 rich media for Mtb, with the
  addition of ADC or OADC growth supplements, and a glycerol carbon source). The
  model files for the iEK1011 and iEK1011_v2 are named as
  iEK1011[\_v2]\_\<media>.json. The simulation models were created for this
  work.
- results: The results from the analysis
  - mtb_tf_results.xlsx: The combined results for the analysis of the Mtb TF
    targets
  - simulation_model_results.xlsx: The combined results for the simulation model
  - simulation: Directory with the individual results for the simulation model
    - iMAT: iMAT results
      - imat_activity.csv: Results of the binary variables for the iMAT
        optimization problem
      - imat_divergence.csv: Divergence (Kullback-Leibler) for all the reactions
        and metabolite synthesis networks between the base simulation model and
        the 'fva' iMAT model
      - escher_maps: Visualizations created with Escher Maps
        - imat_divergence_simulation_model.html: Webpage with the Escher map for
          the divergence between the base simulation model and the 'fva' iMAT
          model
        - imat_solution_simulation_model.html: Webpage with the Escher map
          showing the results of the binary variables in the iMAT optimization
          problem
      - ko_divergence: Results of divergence analysis for all single gene knock
        outs in the simulation model
        - ko_divergence_results.csv: Single gene knock-out divergence
          (Kullback-Leibler) for the subsystems, reactions, and metabolite
          networks in the simulation model
  - mtb_transcription_factors: Directory with the individual results for the Mtb
    TF analysis
    - escher_maps: Directory of various Escher Maps visualizing results
      regarding ArgR, for different parts of the Metabolism
      - ArgR_density\_\<metabolic subsystem>: Escher maps of the gene target
        density of ArgR in different metabolic subsystems
      - ArgR_divergene\_\<metabolic subsystem>: Escher maps of the normalized
        divergence (across TFs) of the ArgR iMAT model
      - ArgR_enrich_pval\_\<metabolic subsystem>: Escher maps of the p-values
        for the neighborood target enrichment of ArgR in various metabolic
        subsystems
    - divergence_results.csv: Divergence of the TFOE iMAT models
    - flux_mi_gene_centrality.csv: Centrality of genes in the flux mutual
      information network (converted from reaction centrality in the network)
    - imat_compare.csv: Analysis of iMAT results for ArgR
    - kegg_divergence.csv: Divergence of KEGG pathways for the TFOE iMAT models
    - kegg_divergence_normalized.csv: The kegg_divergence.csv normalized using a
      StandardScaler from scikit-learn across the various TFs
    - ko_divergence_biomass_statistics.csv: Statistics regarding gene
      essentiality/vulnerability and its relationship to the divergence of the
      BIOMASS\_\_2 reaction following the knock out of the genes
    - ko_divergence_tf_target_analysis.csv: Analysis of the divergence caused by
      the knockout of genes targeted by the TFs
    - metabolic_reaction_network_centrality.csv: Centrality of reactions in the
      stoichiometric connectivity network
    - metabolic_reaction_network_centrality_analysis.csv: Analysis of the
      centrality of reactions targeted by the TFs in the stoichiometric
      connectivity network
    - metabolite_divergence.csv: Divergence of the metabolite sub-networks for
      the TF iMAT models
    - metabolite_divergence_normalized.csv: Divergence of the metabolite
      sub-networks for the TF iMAT models, normalized across all the TF models
    - metabolite_gsva.csv: Results of gene set variation analysis (GSVA) for the
      TFs, using the metabolite networks as the network input for GSVA
    - mutual_information_tf_target_centrality.csv: Results of Mann-Whitney
      U-tests for the centrality, within the flux mutual information network, of
      targets of TFs vs non-targeted genes
    - mutual_information_vi_essentiality_statistics.csv: Statistics comparing
      the centrality of genes, within the flux mutual information network,
      between essential and non-essential genes. Also, includes correlations
      between the mutual information centrality of genes and their Vulnerability
      Index.
    - subsystem_divergence.csv: Divergence of the subsystems in the iEK1011_v2
      model for all of the TFOE iMAT models
    - subsystem_divergence_normalized.csv: Divergence of the subsystems in the
      iEK1011_v2 model for all of the TFOE iMAT models, normalized across TFs
    - tf_target_density.csv: The gene target density of TFs in the
      stoichiometric connectivity network, with rows representing different
      reactions, and columns representing different TFs (and also reaction
      information on the right end of the file)
    - tf_target_enrichment.csv: The neighborhood gene target enrichment of TFs
      in the stoichiometric connectivity network, with rows representing
      different reactions, and columns representing different TFs (and also
      reaction information on the right end of the file)
    - tf_target_metabolite_network_enrichment.csv: Enrichment of the TF gene
      regulatory targets in the metabolite networks of iEK1011_v2.
    - tf_target_subsystem_enrichment.csv: Enrichment of the TF gene regulatory
      targets in the subsystems of iEK1011_v2.
- scripts: Python scripts for performing the analysis
  - collate_results.py: Script for combining the results of the simulation and
    Mtb TF analysis into the final excel files
  - simulation: Directory of scripts for example analysis using a small
    simulated model
    - 01_create_model.py: Script which creates the simulation model
    - 02_metabolite_networks.py: Script for finding the metabolite networks in
      the simulation model
    - 03_metabolic_network_analysis.py: Script for analyzing the stoichiometric
      connectivity network of the simulation model
    - 04_ko_divergence.py: Script for finding the divergence caused by all
      single gene knockouts in the simulation model
    - 05_mutual_information.py: Script for finding the flux mutual information
      network of the simulation model, and evaluating centrality of nodes in
      this network
    - 06_density.py: Script to find the gene target density and enrichment in
      the simulation model
    - 07_imat_simulation.py: Script to create an iMAT model from the simulation
      model and analyze the divergence between this iMAT model and the base
      simulation model
    - 08_visualization.py: Script to generate visualizations for results on the
      simulation model
  - mtb_transcription_factors: Directory of scripts for investigating Mtb TF
    regulation using MetworkPy
    - 00_metabolite_info.py: Script to gather information about the metabolites
      in the iEK1011_v2 model into a single csv file
    - 00_reaction_info.py: Script to gather information about the reactions in
      the iEK1011_v2 model into a single csv file
    - 01_model_generation.py: Script to generate iMAT models for all of the TFOE
      strains
    - 02_model_sampling.py: Script to sample from the TFOE iMAT models
    - 03_divergence.py: Script to calculate the divergence between the TFOE iMAT
      models and the base iEK1011_v2 models
    - 04_tf_ko_divergence.py: Script to find the divergence caused by each
      single gene knockout in the iEK1011_v2 model, and evaluate the gene
      regulatory targets of the TFs for the metabolite subsystems they cause the
      greatest divergence in following knockout
    - 05_tf_target_density.py: Script to find the gene target density and
      enrichment of the TFs in the stoichiometric connectivity network of the
      iEK1011_v2 model
    - 06_reaction_network_centrality.py: Analysis of the centrality of
      reactions within the metabolic stoichiometric connectivity network,
      and the centrality of the TF targets
    - 07_metabolite_gsva.py: Script to perform gene set variation analysis
      (GSVA) on the TF expression data with the metabolite networks as the
      network input of GSVA
    - 08_metabolite_network_target_enrichment.py: Script to analyze the TF
      gene regulatory target enrichment in the metabolite networks of iEK1011_v2.
    - 09_mutual_information_centrality.py: Script to find the centrality of reactions
      in the flux mutual information network, and evaluate the centrality of the
      gene regulatory targets of the TFs
    - 10_mtb_tf_visualization.py: Script for creating visualizations of the results
      of the Mtb TF analyses
    - 11_imat_test.py: Script to perform analysis of ArgR using iMAT
    - common_functions.py: Script with some functions used by the other scripts
      in this directory
</details>

## References

### Software

- COBRApy: Ebrahim, A., Lerman, J.A., Palsson, B.O. et al. COBRApy:
  COnstraints-Based Reconstruction and Analysis for Python. BMC Syst Biol 7, 74
  (2013). [DOI: 10.1186/1752-0509-7-74](https://doi.org/10.1186/1752-0509-7-74)
- Escher Maps: Zachary A. King, Andreas Dräger, Ali Ebrahim, Nikolaus
  Sonnenschein, Nathan E. Lewis, and Bernhard O. Palsson (2015) Escher: A web
  application for building, sharing, and embedding data-rich visualizations of
  biological pathways, PLOS Computational Biology 11(8): e1004321.
  [doi:10.1371/journal.pcbi.1004321](http://dx.doi.org/10.1371/journal.pcbi.1004321)
- NetworkX: Aric A. Hagberg, Daniel A. Schult and Pieter J. Swart, “Exploring
  network structure, dynamics, and function using NetworkX”, in Proceedings of
  the 7th Python in Science Conference (SciPy2008), Gäel Varoquaux, Travis
  Vaught, and Jarrod Millman (Eds), (Pasadena, CA USA), pp. 11–15, Aug 2008
- SciPy: Pauli Virtanen, Ralf Gommers, Travis E. Oliphant, Matt Haberland, Tyler
  Reddy, David Cournapeau, Evgeni Burovski, Pearu Peterson, Warren Weckesser,
  Jonathan Bright, Stéfan J. van der Walt, Matthew Brett, Joshua Wilson, K.
  Jarrod Millman, Nikolay Mayorov, Andrew R. J. Nelson, Eric Jones, Robert Kern,
  Eric Larson, CJ Carey, İlhan Polat, Yu Feng, Eric W. Moore, Jake VanderPlas,
  Denis Laxalde, Josef Perktold, Robert Cimrman, Ian Henriksen, E.A. Quintero,
  Charles R Harris, Anne M. Archibald, Antônio H. Ribeiro, Fabian Pedregosa,
  Paul van Mulbregt, and SciPy 1.0 Contributors. (2020) SciPy 1.0: Fundamental
  Algorithms for Scientific Computing in Python. Nature Methods, 17(3), 261-272.
  DOI: [10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2).

### Data References:

- tfoe_targets.xlsx: Data from 'Rustad, T.R., Minch, K.J., Ma, S. et al. Mapping
  and manipulating the Mycobacterium tuberculosis transcriptome using a
  transcription factor overexpression-derived regulatory network. Genome Biol
  15, 502 (2014).
  [https://doi.org/10.1186/s13059-014-0502-3](https://doi.org/10.1186/s13059-014-0502-3)'
  [PMID: 25380655](https://pmc.ncbi.nlm.nih.gov/articles/PMC4249609/)
- Mycobacterium_tuberculosis_H37Rv_txt_v5.txt: MycoBrowser release V5 from
  'Kapopoulou A, Lew JM, Cole ST. The MycoBrowser portal: a comprehensive and
  manually annotated resource for mycobacterial genomes. Tuberculosis (Edinb).
  Jan 91(1):8-13 (2011)' at
  [https://mycobrowser.epfl.ch/](https://mycobrowser.epfl.ch/)
- bigg_models_metabolites.txt: Information about metabolites in the BiGG
  database, from 'King ZA, Lu JS, Dräger A, Miller PC, Federowicz S, Lerman JA,
  Ebrahim A, Palsson BO, and Lewis NE. BiGG Models: A platform for integrating,
  standardizing, and sharing genome-scale models (2016) Nucleic Acids Research
  44(D1):D515-D522. doi:10.1093/nar/gkv1049' at
  [http://bigg.ucsd.edu/](http://bigg.ucsd.edu/)
- bosch_vi.xlsx: Data from 'Bosch B, DeJesus MA, Poulton NC, Zhang W, Engelhart
  CA, Zaveri A, Lavalette S, Ruecker N, Trujillo C, Wallach JB, Li S, Ehrt S,
  Chait BT, Schnappinger D, Rock JM. Genome-wide gene expression tuning reveals
  diverse vulnerabilities of M. tuberculosis. Cell. 2021 Aug
  19;184(17):4579-4592.e24. doi: 10.1016/j.cell.2021.06.033. Epub 2021 Jul 22.
  PMID: 34297925; PMCID: PMC8382161.'
