<p align="center">
  <img src="Voorblad/Reuma_voorblad.png" alt="Voorblad" width="600"/>
</p>

Transcriptomics onthult immuungerelateerde genexpressie bij Reumatoïde Artritis (RA)

---

Reumatoïde artritis (RA) is een chronische auto-immuunziekte waarbij het immuunsysteem het eigen gewrichtsslijmvlies aanvalt. Dit leidt tot ontsteking, zwelling en uiteindelijk tot beschadiging van kraakbeen en bot. De moleculaire mechanismen die ten grondslag liggen aan de ontwikkeling van RA zijn nog niet volledig opgehelderd, maar eerdere studies wijzen op een centrale rol van immuungerelateerde genen en signaalroutes (zie [RA_gene_expression.pdf](Bronnen/RA_gene_expression.pdf)). 

Eerdere studies (zie [....pdf](Bronnen/.....pdf)) tonen aan dat ontstekingsgerelateerde genen zoals *IL6*, *IL1B* en *TLR4* een sleutelrol spelen bij RA. Ook signaleringsroutes zoals TNF- en Toll-like-receptor pathways zijn herhaaldelijk in verband gebracht met de pathogenese van RA.

RNA-sequencing (RNA-seq) biedt een krachtig middel om de expressie van duizenden genen gelijktijdig te meten. Door het vergelijken van genexpressieprofielen tussen RA-patiënten en gezonde controles kunnen differentieel tot expressie komende genen (DEGs) worden geïdentificeerd. Deze DEGs kunnen vervolgens worden gekoppeld aan biologische processen via Gene Ontology (GO) en aan signaalroutes via KEGG-pathways, wat kan leiden tot nieuwe inzichten in de pathogenese van RA (zie [RNAsequencing.pdf](Bronnen/RNAsequencing.pdf)).

In dit project analyseren we RNA-seq data van synoviaal weefsel afkomstig van vier RA-patiënten en vier gezonde individuen. We identificeren genen die significant op- of neerwaarts gereguleerd zijn, en voeren pathway- en GO-analyse uit om te onderzoeken welke biologische mechanismen mogelijk een rol spelen bij RA. Het doel is om nieuwe aanknopingspunten te vinden voor het begrijpen van de ziekte op transcriptomisch niveau.

--- 

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

## 🔬 Methode

Deze analyse gebruikt RNA-seq data van vier RA-patiënten en vier gezonde controles. De ruwe reads (FASTQ-bestanden) werden uitgelijnd op het humane referentiegenoom (GRCh38) met het `Rsubread`-pakket, waarna `.BAM`-bestanden werden gegenereerd. Deze zijn gesorteerd en geïndexeerd (`Data/processed/`).

Met `featureCounts()` werd een gen-telling uitgevoerd op basis van een GTF-bestand, resulterend in een count-matrix. De differentiële expressie-analyse werd uitgevoerd in `DESeq2`, waarbij log2 fold changes en aangepaste p-waardes (padj) werden berekend. De significante genen (padj < 0.05, |log2FC| > 1) zijn gevisualiseerd in een volcano plot.
Voor functionele interpretatie is een KEGG-pathwayanalyse uitgevoerd met `pathview`. De GO-enrichmentanalyse is uitgevoerd met `goseq`, met biascorrectie via een Probability Weighting Function. De top GO-termen zijn weergegeven in een dot plot.

Alle gebruikte scripts zijn te vinden in [`Scripts/Eigen data.R`](Scripts/Eigen_data.R). Zie het [flowschema](Resultaten/Flowschema.png) voor een overzicht van de workflow.


---


## 📊 Resultaten

Na kwaliteitscontrole en mapping zijn de reads succesvol uitgelijnd op het humane referentiegenoom GRCh38. Op basis van de gen-tellingen (`count_matrix_groot.csv`) werd met `DESeq2` een differentiële expressieanalyse uitgevoerd. Hieruit kwamen meerdere genen significant verschillend tot expressie tussen RA- en controlemonsters.

In de volcano plot ([VolcanoplotWC.png](Resultaten/VolcanoplotWC.png)) zijn deze DEGs visueel weergegeven. Genen met een p-waarde < 0.05 en |log2 fold change| > 1 zijn rood gekleurd; enkele opvallende genen met hoge expressieverandering zijn gelabeld.

Voor GO-analyse is eerst gecorrigeerd voor genlengtebias met een Probability Weighting Function ([pwf_plot.png](Resultaten/pwf_plot.png)). De daadwerkelijke GO-enrichment ([GO_resultaten_plot.png](Resultaten/GO_resultaten_plot.png)) toont dat termen gerelateerd aan immuunrespons, RNA-polymerase II-activiteit sterk verrijkt zijn onder de DE-genen.

De KEGG-pathwayanalyse ([hsa05323.pathview.png](Resultaten/hsa05323.pathview.png)) geeft inzicht in RA-gerelateerde signaalroutes. Genen zoals **IL6**, **IL1B**, en **TLR2/4** vertonen duidelijke opregulatie binnen het ‘Rheumatoid arthritis’ pathway.

Alle resultaten zijn opgeslagen in `Resultaten/`, inclusief de ruwe DESeq2-output (`Resultaten_RA_vs_Normal.csv`) en de gebruikte tellingen (`count_matrix_groot.csv`).




## ✅ Conclusie 

Deze transcriptomics-analyse heeft geleid tot de identificatie van meerdere differentieel tot expressie komende genen (DEGs) tussen RA-patiënten en gezonde individuen. Opvallend waren onder andere genen betrokken bij ontstekingsreacties en immuunactivatie, zoals IL6, TLR2/4 en CXCL-familieleden. De GO-enrichmentanalyse benadrukte vooral termen gerelateerd aan immuunrespons, RNA-polymeraseactiviteit en T-helpercel differentiatie. Deze resultaten zijn in lijn met bestaande literatuur over de rol van het immuunsysteem bij de pathogenese van RA.

De KEGG-pathwayanalyse van het ‘Rheumatoid arthritis’ signaalnetwerk bevestigde deze bevindingen visueel en toonde opregulatie van meerdere ontstekingsgerelateerde routes. De combinatie van genexpressieanalyse en functionele annotatie onderstreept het belang van transcriptomics bij het beter begrijpen van RA op moleculair niveau.

Hoewel het aantal monsters beperkt was (n = 8), bieden de resultaten een waardevolle eerste indruk van genregulatie bij RA. Toekomstig onderzoek zou kunnen uitbreiden met een groter cohort, differentiatie tussen vroege en late RA, en aanvullende celtype-specifieke analyse (bijv. single-cell RNA-seq). Verder validatie via qPCR of proteomics wordt aanbevolen om de biologische relevantie van deze genen in context van RA te bevestigen.

---




## 💻 GitHub en Reproduceerbaarheid

Voor dit project is GitHub ingezet om de RNA-seq analyse reproduceerbaar en transparant te beheren. De repository is logisch gestructureerd met aparte mappen voor ruwe data (data/), scripts (scripts/), resultaten (results/) en documentatie (docs/).

Alle gebruikte scripts zijn stap voor stap gedocumenteerd, inclusief uitleg over welke packages en parameters zijn toegepast. Elke belangrijke analyse (mapping, gen-telling, DE-analyse, KEGG/GO-analyse) is terug te vinden in het script RA_analysis_script.R.

De workflow is volledig reproduceerbaar: met de ruwe FASTQ-bestanden en het script kunnen alle resultaten exact opnieuw worden gegenereerd.

Visualisaties (zoals de volcano plot en KEGG-pathwayplaatjes) zijn opgeslagen als PNG-bestanden en direct toegankelijk via de results/ map.

Versiebeheer via GitHub zorgt ervoor dat iedere aanpassing in het script of data-analyseproces wordt bijgehouden. Zo is transparant te volgen wat er wanneer is aangepast, wat bijdraagt aan de wetenschappelijke betrouwbaarheid van dit project.




