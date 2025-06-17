<p align="center">
  <img src="Voorblad/Reuma_voorblad.png" alt="Voorblad" width="600"/>
</p>

# Transcriptomics Analyse bij Reumatoïde Artritis (RA)

Reumatoïde artritis (RA) is een chronische auto-immuunziekte waarbij het immuunsysteem het eigen gewrichtsslijmvlies aanvalt. Dit leidt tot ontsteking, zwelling en uiteindelijk tot beschadiging van kraakbeen en bot. De moleculaire mechanismen die ten grondslag liggen aan de ontwikkeling van RA zijn nog niet volledig opgehelderd, maar eerdere studies wijzen op een centrale rol van immuungerelateerde genen en signaalroutes [1]

RNA-sequencing (RNA-seq) biedt een krachtig middel om de expressie van duizenden genen gelijktijdig te meten. Door het vergelijken van genexpressieprofielen tussen RA-patiënten en gezonde controles kunnen differentieel tot expressie komende genen (DEGs) worden geïdentificeerd. Deze DEGs kunnen vervolgens worden gekoppeld aan biologische processen via Gene Ontology (GO) en aan signaalroutes via KEGG-pathways, wat kan leiden tot nieuwe inzichten in de pathogenese van RA [2].

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

## 📄 Methode

De analyse begon met het downloaden en uitpakken van RNA-seq data van vier RA-patiënten en vier gezonde controles. De reads werden uitgelijnd op het humane referentiegenoom (GRCh38) met het Rsubread-pakket in R, waarna de resulterende .BAM-bestanden werden gesorteerd en geïndexeerd. Vervolgens werd met featureCounts() een telling per gen uitgevoerd op basis van een GTF-annotatiebestand.

De gegenereerde count-matrix werd ingelezen in het DESeq2-pakket om differentieel tot expressie komende genen (DEGs) te identificeren. Hierbij werd per gen de log2 fold change en aangepaste p-waarde berekend. Significante DEGs (padj < 0.05 en |log2FC| > 1) werden gevisualiseerd in een volcano plot.

Voor functionele interpretatie zijn KEGG-pathways geanalyseerd met het pathview-pakket. Gen-ID’s werden gemapt op het ‘Rheumatoid Arthritis’ pathway (hsa05323). Daarnaast is een GO-enrichmentanalyse uitgevoerd met het goseq-pakket, waarbij werd gecorrigeerd voor genlengtebias (PWF). De top GO-termen zijn gevisualiseerd in een bubbleplot.

▸ De gebruikte scripts zijn te vinden in [`Scripts/`](Scripts/)  
▸ De gegenereerde resultaten (tabellen en figuren) staan in [`Resultaten/`](Resultaten/)  
▸ Alle verwerkte inputdata (zoals `.BAM`-bestanden, tellingen en annotaties) zijn beschikbaar in [`Data/processed/`](Data/processed/)   
▸ Zie het [flowschema](Resultaten/method_flowchart.png) voor een visueel overzicht van de analysepipeline 

---


## 📊 Resultaten




## ✅ Conclusie 

Op basis van de RNA-seq analyse zijn meerdere genen gevonden die significant verschillen in expressie tussen RA-patiënten en gezonde controlepersonen. Deze genen zijn vooral betrokken bij immuungerelateerde processen, zoals cytokinesignalering, ontstekingsreacties en leukocytenactivatie.

De KEGG- en GO-analyse bevestigen dat reumatoïde artritis gepaard gaat met verhoogde activiteit van ontstekingsroutes. Deze bevindingen ondersteunen het beeld dat RA een chronische ontstekingsziekte is met een sterke transcriptomische footprint in synoviaal weefsel.

Deze studie laat zien dat RNA-seq een krachtig hulpmiddel is om inzicht te krijgen in de moleculaire mechanismen van RA en biedt mogelijke aanknopingspunten voor biomarkerontwikkeling of gerichte therapieën.

### 📚 Bronnen

[1] [RA_gene_expression.pdf](Bronnen/RA_gene_expression.pdf) – overzicht van genexpressie bij RA    
[2] [RNA_sequencing.pdf](Bronnen/RNAsequencing.pdf) – achtergrond over pathway- en GO-analyse


## 💻 GitHub en Reproduceerbaarheid

Voor dit project is GitHub ingezet om de RNA-seq analyse reproduceerbaar en transparant te beheren. De repository is logisch gestructureerd met aparte mappen voor ruwe data (data/), scripts (scripts/), resultaten (results/) en documentatie (docs/).

Alle gebruikte scripts zijn stap voor stap gedocumenteerd, inclusief uitleg over welke packages en parameters zijn toegepast. Elke belangrijke analyse (mapping, gen-telling, DE-analyse, KEGG/GO-analyse) is terug te vinden in het script RA_analysis_script.R.

De workflow is volledig reproduceerbaar: met de ruwe FASTQ-bestanden en het script kunnen alle resultaten exact opnieuw worden gegenereerd.

Visualisaties (zoals de volcano plot en KEGG-pathwayplaatjes) zijn opgeslagen als PNG-bestanden en direct toegankelijk via de results/ map.

Versiebeheer via GitHub zorgt ervoor dat iedere aanpassing in het script of data-analyseproces wordt bijgehouden. Zo is transparant te volgen wat er wanneer is aangepast, wat bijdraagt aan de wetenschappelijke betrouwbaarheid van dit project.




