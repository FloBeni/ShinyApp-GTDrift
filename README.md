[![DOI](https://zenodo.org/badge/706581407.svg)](https://zenodo.org/doi/10.5281/zenodo.10022520)

This archive includes a Shiny Application to explore the database from the paper 'GTDrift: A resource for exploring the interplay between genetic drift, genomic and transcriptomic characteristics in eukaryotes'.

Developed by Florian Bénitière, Laurent Duret and Anamaria Necsulea. Laboratoire de Biométrie et Biologie Évolutive, Université Lyon 1, UMR CNRS 5558, Villeurbanne, France.

A light version of the database is present in 'www/database' which corresponds to database downloadable at [![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.10908656.svg)](https://doi.org/10.5281/zenodo.10908656).

The application is available at: <https://lbbe-shiny.univ-lyon1.fr/ShinyApp-GTDrift/>

There are 6 tabs.

-   The tab 'Inter-species graphics' facilitates the comparison of genome characteristics across different species through graphical representation. Additionally, users have the option to upload data in a tab-separated text format, where each species is represented in a separate row, with the variables of interest organized in columns. An illustrative example of such a tabular dataset can be found in the table accessible at 'www/species_informations_tables/data_by_species.tab'. Important: the provided table must have a column 'species' with the scientific name of the species (*i.e.* Arabidopsis_thaliana).

-   The tab 'Inter-species Axis' described the variables that are available in the 'Inter-species graphics' tab.

-   The tab 'Intra-species graphics' permit to explore characteristics within a species between introns or between genes. Also, users can download metadata for BUSCO annotation, genes expression profile or introns splicing events (see method in original paper).

-   The tab 'Intra-species Axis' described the variables that are proposed in the 'Intra-species graphics' tab.

-   In the 'Gene structure' tab, you can explore the introns detected in RNA-seq alignments for a specific gene. The introns are color-coded based on categories such as their location within the CDS or outside of it, as well as whether they are major-isoform or minor-isoform introns (as described in the original paper's methodology).

-   'Phylogenetic tree' tab facilitates the examination of phylogenetic trees used for Phylogenetic Generalized Least Squares regression within 'Inter-species graphics' tab.
