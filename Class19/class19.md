# Class19
Mason Lew (PID: A17533139)

> Q1. \[1pt\]What protein do these sequences correspond to? (Give both
> full gene/protein name and official symbol).

``` r
library(bio3d)
data <- read.fasta("/Users/Mason/Desktop/BIMM143/BIMM-142/Class19/A17533139_mutant_seq.fa")
data
```

                   1        .         .         .         .         .         60 
    wt_healthy     MGMACLTMTEMEGTSTSSIYQNGDISGNANSMKQIDPVLQVYLYHSLGKSEADYLTFPSG
    mutant_tumor   MGMACLTMTEMEGTSTSSIYQNGDISGNANSMKQIDPVLQVYLYHSLGKSEADYLTFPSG
                   ************************************************************ 
                   1        .         .         .         .         .         60 

                  61        .         .         .         .         .         120 
    wt_healthy     EYVAEEICIAASKACGITPVYHNMFALMSETERIWYPPNHVFHIDESTRHNVLYRIRFYF
    mutant_tumor   EYVAEEICIAASKACGITPVYHNMFALMSETERIWYPPNHVFHIDESTRHNVLYRIRFYF
                   ************************************************************ 
                  61        .         .         .         .         .         120 

                 121        .         .         .         .         .         180 
    wt_healthy     PRWYCSGSNRAYRHGISRGAEAPLLDDFVMSYLFAQWRHDFVHGWIKVPVTHETQEECLG
    mutant_tumor   PRWYCSGSNRAYRHGISRGAEAPLLDDFVMSYLFAQWRHDFVHGWIKVPVTHETQEECLG
                   ************************************************************ 
                 121        .         .         .         .         .         180 

                 181        .         .         .         .         .         240 
    wt_healthy     MAVLDMMRIAKENDQTPLAIYNSISYKTFLPKCIRAKIQDYHILTRKRIRYRFRRFIQQF
    mutant_tumor   MAVLDMMRIAKENDQTPLAIYNSISYKTFLPKCIRAKIQDYHILTRKRIRYRFRRFIQQF
                   ************************************************************ 
                 181        .         .         .         .         .         240 

                 241        .         .         .         .         .         300 
    wt_healthy     SQCKATARNLKLKYLINLETLQSAFYTEKFEVKEPGSGPSGEEIFATIIITGNGGIQWSR
    mutant_tumor   SQCKATARNLKLKYLINLETLQSAFYTEKFEVKEPGSGPSGEEIFATIIITGNGGIQWSR
                   ************************************************************ 
                 241        .         .         .         .         .         300 

                 301        .         .         .         .         .         360 
    wt_healthy     GKHKESETLTEQDLQLYCDFPNIIDVSIKQANQEGSNESRVVTIHKQDGKNLEIELSSLR
    mutant_tumor   GKHKESETLTEQDLQLYCDFPNIIDVSIKQANQEGSNESRVVTIHKQDGKNLEIELSSLR
                   ************************************************************ 
                 301        .         .         .         .         .         360 

                 361        .         .         .         .         .         420 
    wt_healthy     EALSFVSLIDGYYRLTADAHHYLCKEVAPPAVLENIQSNCHGPISMDFAISKLKKAGNQT
    mutant_tumor   EALSFVSLIDGYYRLTADAHHYLCKEVAPPAVLENIQSNCHGPISMDFAISKLKKAGNQT
                   ************************************************************ 
                 361        .         .         .         .         .         420 

                 421        .         .         .         .         .         480 
    wt_healthy     GLYVLRCSPKDFNKYFLTFAVERENVIEYKHCLITKNENEEYNLSGTKKNFSSLKDLLNC
    mutant_tumor   GLYVLRCSPKDFNKYFLTFAVERENVIEYKHCLITKNENEEYNLSGTKKNFSSLKDLLNC
                   ************************************************************ 
                 421        .         .         .         .         .         480 

                 481        .         .         .         .         .         540 
    wt_healthy     YQMETVRSDNIIFQFTKCCPPKPKDKSNLLVFRTNGVSDVPTSPTLQRPTHMNQMVFHKI
    mutant_tumor   YQMETVRSDNIIFQFTKCCPPKPKDKSNLLVFRTNGVSDVPTSPTLQRPTHMNQMVFHKI
                   ************************************************************ 
                 481        .         .         .         .         .         540 

                 541        .         .         .         .         .         600 
    wt_healthy     RNEDLIFNESLGQGTFTKIFKGVRREVGDYGQLHETEVLLKVLDKAHRNYSESFFEAASM
    mutant_tumor   RNEDLIFNESLGQGTFVKIFKGVRREVGDYGQLHETEVLLKVLDKAHRNYEESFFEAASM
                   **************** ********************************* ********* 
                 541        .         .         .         .         .         600 

                 601        .         .         .         .         .         660 
    wt_healthy     MSKLSHKHLVLNYGVCVCGDENILVQEFVKFGSLDTYLKKNKNCINILWKLEVAKQLAWA
    mutant_tumor   MSRLSHKHLVLNYGVCVCGDENILVQEFVKFGSLDTYLKKNKNCINILWKLEVAKQLAWA
                   **^********************************************************* 
                 601        .         .         .         .         .         660 

                 661        .         .         .         .         .         720 
    wt_healthy     MHFLEENTLIHGNVCAKNILLIREEDRKTGNPPFIKLSDPGISITVLPKDILQERIPWVP
    mutant_tumor   MHFLEENTLIHGNVCAKNILLIREEDRKYGNPPFIKLSDPGISITVLPKDILQERIPWVP
                   **************************** ******************************* 
                 661        .         .         .         .         .         720 

                 721        .         .         .         .         .         780 
    wt_healthy     PECIENPKNLNLATDKWSFGTTLWEICSGGDKPLSALDSQRKLQFYEDRHQLPAPKWAEL
    mutant_tumor   PECIENPKNLNLATDKWSFGTTLWEICSGGDKPLSALDSQRKLQFYEDRHQLPAPKWAEL
                   ************************************************************ 
                 721        .         .         .         .         .         780 

                 781        .         .         .         .         .         840 
    wt_healthy     ANLINNCMDYEPDFRPSFRAIIRDLNSLFTPDYELLTENDMLPNMRIGALGFSGAFEDRD
    mutant_tumor   ANLINNCMDYEPDFRPSFRAIIRDLNSLFTPDYELLTENDMLPNMRIGALGFSGAFEDRD
                   ************************************************************ 
                 781        .         .         .         .         .         840 

                 841        .         .         .         .         .         900 
    wt_healthy     PTQFEERHLKFLQQLGKGNFGSVEMCRYDPLQDNTGEVVAVKKLQHSTEEHLRDFEREIE
    mutant_tumor   PTQFEERHLKFLQQLGKGNFGSVEMCRYDPLQDNTGEVVAVKKLQHSTEEHLRDFEREIE
                   ************************************************************ 
                 841        .         .         .         .         .         900 

                 901        .         .         .         .         .         960 
    wt_healthy     ILKSLQHDNIVKYKGVCYSAGRRNLKLIMEYLPYGSLRDYLQKHKERIDHIKLLQYTSQI
    mutant_tumor   ILKSLQHDNIVKYKGVCYSAGRRNLKLIMEYLPYGSLRDYLQKHKERIDHIKLLQYTSQI
                   ************************************************************ 
                 901        .         .         .         .         .         960 

                 961        .         .         .         .         .         1020 
    wt_healthy     CKGMEYLGTKRYIHRDLATRNILVENENRVKIGDFGLTKVLPQDKEYYKVKEPGESPIFW
    mutant_tumor   CKGMEYLGTKRYIHRDLATRNILVENENRVKIGDFGLTKVLPQDKEYYKVKEPGESPIFW
                   ************************************************************ 
                 961        .         .         .         .         .         1020 

                1021        .         .         .         .         .         1080 
    wt_healthy     YAPESLTESKFSVASDVWSFGVVLYELFTYIEKSKSPPAEFMRMIGNDKQGQMIVFHLIE
    mutant_tumor   YAPESLTESKFSVASDVWSFGVVLYELFTYIEKSKSPPAEFMRMIGNDKQGQMIVFHLIE
                   ************************************************************ 
                1021        .         .         .         .         .         1080 

                1081        .         .         .         .         . 1132 
    wt_healthy     LLKNNGRLPRPDGCPDEIYMIMTECWNNNVNQRPSFRDLALRVDQIRDNMAG
    mutant_tumor   LLKNNGRLPRPDGCPDEIYMIMTECWNNNVNQRPSFRDLALRVDQIRDNMAG
                   **************************************************** 
                1081        .         .         .         .         . 1132 

    Call:
      read.fasta(file = "/Users/Mason/Desktop/BIMM143/BIMM-142/Class19/A17533139_mutant_seq.fa")

    Class:
      fasta

    Alignment dimensions:
      2 sequence rows; 1132 position columns (1132 non-gap, 0 gap) 

    + attr: id, ali, call

tyrosine-protein kinase JAK2 isoform a \[Homo sapiens\]

``` r
conserve <- conserv(data)
positions <- which(conserve != 1)
positions
```

    [1] 557 591 603 689

> Q2.
> \[6pts\]Whatarethetumorspecificmutationsinthisparticularcase(e.g.A130V)?

``` r
residues <- data$ali[positions]
residues
```

    [1] "P" "I" "K" "H"

T557V, S591E, K603R, T689Y

> Do your mutations cluster to any particular domainand if so give the
> name and PFAM id of this domain? Alternately note whether your protein
> is single domainand provide it’s PFAM id/accession and name
> (e.g.PF00613 and PI3Ka).

This protein is not a single domain, the mutations are focused within
the Protein tyrosine and serine/threonine kinase domain: PFAM id:
PF07714, Accession: A0A010QA77

> Q4. \[2pts\]Using the NCI-GDC list the observed top 2 missense
> mutations in this protein (aminoacidsubstitutions)?

R683G, T875N

> Q5. \[2pts\]What two TCGA projects have the most cases affected by
> mutations of this gene? (Give the TCGA “code” and “ProjectName” for
> example “TCGA-BRCA”and “Breast Invasive Carcinoma”).

TCGA-UCEC, Uterine Corpus Endometrial Carcinoma

TCGA-STAD, Stomach Adenocarcinoma

> Q6. \[3pts\]List one RCSB PDB identifier with 100% identity to the
> wt_healthy sequence and detail the percent coverage of your query
> sequence for this known structure? Alternately, provide the most
> similar in sequence PDB structure along with it’s percent identity,
> coverage and E-value. Does this structure “cover” (i.e. include or
> span the amino acid residue positions) of your previously identified
> tumor specific mutations?

tyrosine-protein kinase JAK2 isoform a \[Homo sapiens\], NP_001309123.1:
100% query coverage, 100% identity, 0 E-value. Yes this sequence cover
the identified tumor specific mutation

> Q7. \[10pts\] Using Alpha Fold notebook generate a structural model
> using the default parameters for your mutant sequence.

![“/Users/Mason/Desktop/BIMM143/BIMM-142/Class19/587E42C122BFB93B.png”](hi.png)

``` r
cif.file<- "/Users/Mason/Desktop/BIMM143/BIMM-142/Class19/fold_mutant_jak2_model_4.cif"
x <- read.cif(cif.file)
write.pdb(x, file="FT_map5.pdb")
```
