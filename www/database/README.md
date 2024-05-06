This archive include GTDrift database described in the study 'GTDrift: A resource for exploring the interplay between genetic drift, genomic and transcriptomic characteristics in eukaryotes'. The database is generated from [![DOI](https://zenodo.org/badge/685027688.svg)](https://zenodo.org/doi/10.5281/zenodo.10022493).

Developed by Florian Bénitière, Laurent Duret and Anamaria Necsulea. Université de Lyon, Université Lyon 1, CNRS, Laboratoire de Biométrie et Biologie Évolutive UMR 5558, F-69622 Villeurbanne, France.

### Content of the database

-   The file **'list_species.tab'** contains for each species of the database.

    <div>

    -   **species**: Scientific name of the species.

    -   **NCBI.taxid**: taxID from the NCBI.

    -   **assembly_accession**: Assembly accession used.

    -   **clade**: Clade corresponding to the species, it can be the class for vertebrates or the order for insects.

    -   **clade_group**: Clade corresponding to the species, mostly the rank.

    -   **expression_data**: Does the species have transcriptomic data.

    -   **dnds_data**: Does the species have molecular rate (dN/dS) data.

    -   **lht_data**: Does the species have life history traits.

    -   **Ne_data**: Does the species have a polymorphism-derived Ne estimate.

    -   **nb_rnaseq**: Number of RNA-seq studied.

    -   **entire_taxonomy**: Concatenation of the NCBI taxonomy, collected from the R package taxize.

    </div>

-   The file **'life_history_traits_and_polymorphism_derived_Ne.tab'** contains the maximum values found for each species and each life history traits, along with polymorphism-derived Ne estimates.

    <div>

    -   **species**: Scientific name of the species.

    -   **NCBI.taxid**: taxID from the NCBI.

    -   **mass_kg**: Maximum body mass value found (kg).

    -   **mass_db**: Where the body mass value was found (Web database, references...).

    -   **lifespan_days**: Maximum lifespan value found (days).

    -   **lifespan_db**: Where the lifespan value was found (Web database, references...).

    -   **length_cm**: Maximum body length value found (cm).

    -   **length_db**: Where the body length value was found (Web database, references...).

    -   **polymorphism_derived_Ne**: Polymorphism-derived Ne estimates.

    </div>

-   The **'Transcriptomic'** directory contains processed data in tab-separated text format, which corresponds to splicing analyses.

    -   Each species possesses a directory '*`Scientific name`*\_NCBI.taxid*`taxID`*', within can be found a subdirectory with the *`assembly accession`*.

        -   'by_gene_analysis.tab.gz' which contain per-gene data.

        <div>

        -   The header of the GFF used followed by a brief description of the main variables.

        -   **gene_id**: ID of the gene found in the GFF.

        -   **gene_name**: Corresponds to the attribute 'gene=' in the GFF.

        -   **seq_id**: Refers to seqname in the GFF, corresponding to the name of the chromosome or scaffold.

        -   **start**: Start position.

        -   **end**: End position.

        -   **strand**: Defined as + (forward) or - (reverse).

        -   **type**: Feature type name, e.g. Gene, Variation, Similarity found in the GFF.

        -   **attributes**: Semicolon-separated list of tag-value pairs, providing additional information about each feature, extracted from the GFF.

        -   **weighted_fpkm**: Computed as the average FPKM across samples, weighted by the sequencing depth of each sample. The sequencing depth of a sample is the median *per*-base read coverage across BUSCO genes.

        -   **median_fpkm**: Median FPKM across all RNA-seq samples.

        -   **mean_fpkm**: Average FPKM across all RNA-seq samples.

        -   **std_fpkm**: Standard deviation FPKM across all RNA-seq samples.

        -   **exon_coverage**: Exonic read *per* bp (measured on all RNA-seq pooled).

        </div>

        -   'by_intron_analysis.tab.gz' which contain per-intron data.

        <div>

        -   **gene_id**: ID of the gene found in the GFF.

        -   **seqname**: seqname in the GFF, corresponding to the name of the chromosome or scaffold

        -   **strand**: Defined as 1 (forward) or -1 (reverse).

        -   **splice5**: Corresponds to the 5' splice donor site.

        -   **splice3**: Corresponds to the 3' splice acceptor site.

        -   **ns**: Number of spliced reads corresponding to the precise excision of the focal intron ($N_s$ in the paper).

        -   **na_spl5**: Number of reads corresponding to alternative splice variants relative to this intron (*i.e.* sharing 5' splice donor site) ($\mathrm{na\_spl5+na\_spl3=N_a}$ in the paper).

        -   **na_spl3**: Number of reads corresponding to alternative splice variants relative to this intron (*i.e.* sharing 3' splice acceptor site) ($\mathrm{na\_spl5+na\_spl3=N_a}$ in the paper).

        -   **nu_spl3**: Number of unspliced reads, co-linear with the genomic sequence (at the 3' splice acceptor site) ($\mathrm{nu\_spl5+nu\_spl3=N_u}$ in the paper).

        -   **nu_spl5**: Number of unspliced reads, co-linear with the genomic sequence (at the 5' splice donor site). ($\mathrm{nu\_spl5+nu\_spl3=N_u}$ in the paper).

        -   **splice_variant_rate**: Proportion of reads alternatively spliced, $\mathrm{1-RAS=AS=\frac{N_a}{N_s~+~N_a}}$.

        -   **nonsplice_variant_rate**: Proportion of unspliced reads, $\mathrm{1-RANS=\frac{N_u}{2\times N_s~+~N_u}}$.

        -   **intron_class**: Three categories of introns: major-isoform introns, defined as those introns that have RANS $>$ 0.5 and RAS $>$ 0.5; minor-isoform introns, defined as those introns that have RANS $\leq$ 0.5 or RAS $\leq$ 0.5; unclassified introns, which do not satisfy the above conditions.

        -   **into_cds**: Check if intron is located within protein-coding regions. To do this, for each protein-coding gene, we extracted the start codons and the stop codons for all annotated isoforms. We then identified the minimum start codon and the maximum end codon positions and we excluded introns that were upstream or downstream of these extreme coordinates.

        -   **id**: Semicolon-separated list of tag-value seqname;gene_id;splice3;end:strand

        -   **fpkm**: Weighted_fpkm of the corresponding gene. Computed as the average FPKM across samples, weighted by the sequencing depth of each sample. The sequencing depth of a sample is the median *per*-base read coverage across BUSCO genes.

        -   **splicesite**: Identify the dinucletotides donor and acceptor (ex: GT AG).

        -   **Annotation**: Check if the intron is annotated (*i.e.* present in the GFF) or not.

        </div>

        -   'SRAruninfo.tab' which is the runinfo extracted from SRA containing information on each RNA-seq Run.

        -   In the directory 'Run' for each Run there is two tables 'by_gene_db.tab.gz' and 'by_intron_db.tab.gz'. These two tables corresponds in line to the previous tables 'by_gene_analysis.tab.gz' and 'by_intron_analysis.tab.gz' respectively.

            -   'by_gene_db.tab.gz' contains the exon_coverage and the expression (fpkm) of the specific Run.

            -   'by_intron_db.tab.gz' contains $N_s$, $N_a$ and $N_u$ variables calculated in the specific Run.

-   The directory labeled **'dNdS'** houses curated data pertaining to the molecular evolutionary rate, specifically, the dN/dS ratio. This ratio signifies the rate at which non-synonymous substitutions occur in comparison to synonymous substitutions. It encompasses four distinct methods for estimation, each associated with specific taxonomic categories: 'Embryophyta', 'Metazoa', 'Eukaryota'. In an effort to encompass a broader spectrum of our species, a fourth methodology was employed named 'per clade', involving the division of clades within the metazoa taxon (see paper for methods).

    -   The file '\*.tab' contains the summarized molecular rate at terminal branches for each species.

        <div>

        -   **species**: Scientific name of the species.

        -   **NCBI.taxid**: taxID from the NCBI.

        -   **dN**: The non-synonymous substitutions rate.

        -   **dS**: The synonymous substitutions rate.

        -   **dNdS**: The rate at which non-synonymous substitutions occur in comparison to synonymous substitutions. A proxy of the effective population size.

        </div>

    -   The directory 'newick' contains for each analysis the dN and dS estimated over the all phylogenetic tree.

    -   The directory 'phylogeny' contains the inferred phylogenetic trees rooted using TimeTree (Kumar 2022).

-   In the provided directory **'BUSCO_annotations'**, we offer the corresponding single-copy orthologous BUSCO gene id for each gene id, derived from the processed BUSCO datasets. The utilized BUSCO datasets encompass three categories: eukaryota, metazoa, and embryophyta.

    -   Each species possesses a directory '*`Scientific name`*\_NCBI.taxid*`taxID`*', within can be found a subdirectory with the *`assembly accession`*.

        -   The file 'busco_to_gene_id\_*`busco set`*.gz' busco_id status sequence score length gene_id duplicated_id cds_length

        <div>

        -   **busco_id**: The BUSCO ID.

        -   **status**: The BUSCO status.

        -   **sequence**: Corresponding protein sequences, we kept the longest for a given gene.

        -   **score**: The BUSCO score.

        -   **length**: Length of the alignment.

        -   **gene_id**: Gene that corresponds to the protein.

        -   **duplicated_id**: BUSCO id followed by gene id, this metric is employed to determine the unambiguous association of a BUSCO gene to a single annotated gene and vice versa.

        -   **cds_length**: Length of the annotated coding sequences (from GFF).

        </div>

