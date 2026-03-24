# MetworkPy Application Note

Data and Code associated with 'MetworkPy: A Python Package for Graph- and
Information-theoretic Investigation of Metabolic Networks'.

## Directory Structure

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
  [A systematic evaluation of Mycobacterium tuberculosis Genome-Scale Metabolic Networks](https://pmc.ncbi.nlm.nih.gov/articles/PMC7316355/).
  The simulation models were created for this work.

## References

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
