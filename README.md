<p align="center">
  <img src="Voorblad/Reuma_voorblad.png" alt="Voorblad" width="600"/>
</p>

# Transcriptomics onthult immuungerelateerde genexpressie bij Reumatoïde Artritis (RA)

## 📁 Inhoud/structuur

- `Data/raw/` – Ruwe datasets die als input dienen voor de analyse. 
- `Data/processed` - Verwerkte datasets die zijn gegenereerd met behulp van scripts.
- `Scripts/` – R-scripts voor preprocessing, analyse en visualisatie.
- `Resultaten/` - Grafieken, tabellen en outputbestanden uit de analyses.
- `Bronnen/` - Externe bronnen zoals artikelen, handleidingen of databronnen. 
- `README.md` - Deze pagina met een toelichting op de projectstructuur en inhoud.
- `Voorblad/` - Documenten die horen bij de opmaak van het verslag
- `data_stewardship/` - Mappen en bestanden om aan te tonen hoe je de projectdata beheert.

---
## Introductie 

Reumatoïde artritis (RA) is een chronische auto-immuunziekte waarbij het immuunsysteem het eigen gewrichtsslijmvlies aanvalt. Dit leidt tot ontsteking, zwelling en uiteindelijk tot beschadiging van kraakbeen en bot. Hoewel de exacte oorzaak onbekend is, is er de afgelopen jaren veel vooruitgang geboekt in de behandeling ([Radu et al., 2021](Bronnen/Radu_2021_RA_management.pdf)).

Dankzij RNA-sequencing zijn er belangrijke genexpressiepatronen in RA ontdekt. Onderzoeken tonen aan dat genen betrokken bij het celskelet, evenals miRNA's, duidelijk anders tot expressie komen bij vroege RA in vergelijking met gezonde personen ([Platzer et al., 2019](Bronnen/Platzer_2019_RA_gene_expression.pdf)). RNA-sequencing maakt het mogelijk om differentieel tot exressie komende genen (DEGs) te indentificeren en deze te koppelen aan biologische processen zoals Gene Ontology (GO) en KEGG-pathways. Deze aanpak geeft meer inzicht in de moleculaire mechanismen van RA ([Zhang et al., 2022](Bronnen/Zhang_2022_RNA_sequencing.pdf)).

In dit project wordt met behulp van transcriptomische analyse onderzocht welke genen bij RA anders tot expressie komen dan bij gezonde personen en welke biologische processen daarbij betrokken zijn. Het doel is om via deze analyse meer inzicht te krijgen in de werking van RA op genetisch niveau.


## 🔬 Methode

Voor deze analyse is RNA-sequence data gebruikt van vier RA-patiënten en vier gezonde controles. De ruwe reads [FASTQ-bestanden](Data/Raw) zijn afkomstig uit een eerder gepubliceerde studie ([Platzer et al., 2019](Bronnen/Platzer_2019_RA_gene_expression.pdf)). Een overzicht van deze samples is te vinden in [Sample_metadata](Data/Raw/sample_metadata_table_RA.png). Deze zijn uitgelijnd op het humane referentiegenoom GRCh38.p14 (NCBI RefSeq: GCF_000001405.40) met behulp van het `Rsubread`-pakket in  [01_preprocessing_alignment](Scripts/01_preprocessing_alignment.R), waarna `.BAM`-bestanden zijn gegenereerd. Deze zijn gersorteerd en geïndexeerd in [02_sort_index_counts](Scripts/02_sort_index_counts.R) en opgelsagen in [Data/Processed](Data/Processed). Daaarnaast is met  met `featureCounts()` een gen-telling uitgevoerd op basis van een GTF-bestand, resulterend in een count-matrix. 

Vervolgens is met `DESeq2` in [03_deseq2_analysis_volcano](Scripts/03_deseq2_analysis_volcano.R) een differentiële expressieanalyse uitgevoerd, waarbij log2 fold changes en aangepaste p-waardes (padj) zijn berekend. De significante genen (padj < 0.05, |log2FC| > 1) zijn gevisualiseerd in een volcano plot.

Voor functionele interpretatie is in [04_kegg_pathway](Scripts/04_kegg_pathway.R) een KEGG-pathwayanalyse uitgevoerd met pathview, gericht op het 'Rheumatoid Arthritis' pathway. In [05_go_enrichment](Scripts/05_go_enrichment.R) is een GO-enrichmentanalyse uitgevoerd met goseq, inclusief biascorrectie via een Probability Weighting Function (pwf). De top GO-termen zijn weergegeven in een dot plot.

Een volledig overzicht van de workflow is te vinden in het [Flowschema](Resultaten/Flowschema.png) . Het samengevoegde script waarin alle stappen zijn opgenomen, is beschikbaar in [Volledige_script](Scripts/Volledige_script.R).

## 📊 Resultaten

Uit de RNA-seq-analyse zijn meerdere genen gevonden die significant verschillen in expressie tussen reumatoïde artritis (RA)-patiënten en gezonde controles.

In de volcano plot ([Volcanoplot](Resultaten/VolcanoplotWC.png)) zijn deze DEGs visueel weergegeven. Genen met een p-waarde < 0.05 en |log2 fold change| > 1 zijn rood gekleurd. Een opvallend sterk gereguleerd gen dat geassocieerd wordt met RA is MT-ND6.

Voor GO-analyse is eerst gecorrigeerd voor genlengtebias met een Probability Weighting Function ([pwf_plot](Resultaten/pwf_plot.png)). De GO-enrichment ([GO_resultaten_plot](Resultaten/GO_resultaten_plot.png)) laat zien dat vooral processen zoals immuunrespons en RNA-polymerase II-activiteit verrijkt zijn bij de genen die verschillend tot expressie komen

De KEGG-pathwayanalyse toonde verhoogde activiteit binnen het ‘Rheumatoid arthritis’ pathway ([hsa05323.pathview](Resultaten/hsa05323.pathview.png)). In het KEGG-diagramGenen zijn meerdere genen betrokken bij ontstekingsroute, waaronder **IL6**, **IL1B** en **TLR2/4**, die duidelijk opgereguleerd zijn (rood).

Zie [Genen_literatuur_tabel](Bronnen/Genen_literatuur_tabel.xlsx)) voor toelichting op de literatuurverwijzingen van specifieke genen.

## ✅ Conclusie 

Deze RNA-seq analyse heeft geleid tot de identificatie van meerdere genen die significant verschillend tot expressie komen tussen RA-patiënten en gezonde controles. Opvallende genen zoals IL6, IL1B, TLR2/4 waren sterk opgereguleerd en zijn bekend om hun rol in ontstekings- en immuunprocessen. De KEGG-pathwayanalyse bevestigde deze bevindingen visueel en toonde verhoogde activiteit binnen het ‘Rheumatoid arthritis’ signaalpad, met duidelijke opregulatie van ontstekingsgerelateerde genen.

Het mitochondriale gen MT-ND6 was significant opgereguleerd. Dit komt overeen met eerder onderzoek waarin mitochondriale peptiden, zoals MT-ND6, betrokken zijn bij neutrofielgemedieerde inflammatie in RA ([Duvvuri et. all, 2021](Bronnen/Duvvuri_2021_MT-ND6.pdf)). 

Hoewel het aantal monsters beperkt was (n = 8), biedt deze analyse een waardevol eerste inzicht in genregulatie bij RA. Vervolgonderzoek met grotere patiëntengroepen en meer gedetailleerde celtypespecifieke methoden zoals single-cell RNA-seq kan bijdragen aan een beter begrip van de betrokken genen en bevestiging van deze bevindingen.





