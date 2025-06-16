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

## ⚙️ Methode
📁 Dataset

De dataset bevat RNA-seq gegevens van synoviumbiopten van 4 gezonde controlepersonen en 4 patiënten met reumatoïde artritis (RA). De ruwe data zijn afkomstig van de Sequence Read Archive (SRA) en bestaan uit paired-end FASTQ-bestanden (~40.000 reads per sample). De sequentiegegevens zijn gegenereerd met een Illumina-platform (exacte platforminformatie: zie metadata bij SRA-toegang).

🔬 Preprocessing & Mapping

De FASTQ-bestanden zijn uitgelijnd (gemapt) tegen het humane referentiegenoom GRCh38 (GCF_000001405.40) met behulp van het align()-commando uit het Rsubread-pakket. Per sample zijn de gegenereerde BAM-bestanden gesorteerd en geïndexeerd met Rsamtools. Voor de genannotatie is gebruik gemaakt van een bijbehorend GTF-bestand (GRCh37.p13), waarbij enkel exon-regio’s werden geteld.

align(index = "ref_human",
      readfile1 = "sample_1.fastq",
      readfile2 = "sample_2.fastq",
      output_file = "sample1.BAM")

📊 Genexpressiematrix

De uitgelijnde reads werden geteld met featureCounts() uit het Rsubread-pakket. Genexpressiewaarden zijn berekend per gen op basis van samengevoegde exon-regio’s. De gegenereerde count-matrix werd opgeslagen als CSV voor verdere analyse.

📈 Differentiële genexpressieanalyse

De genexpressiematrix werd ingelezen in DESeq2. Met behulp van een model waarin conditie (RA of gezond) als designfactor werd opgenomen, is een differentiële expressieanalyse uitgevoerd. De genen werden beschouwd als significant bij een aangepaste p-waarde (FDR) < 0.05 en een absolute log2 fold change > 1.

dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = treatment_table,
                              design = ~ treatment)
dds <- DESeq(dds)
resultaten <- results(dds)

🧬 KEGG pathway-analyse

De differentieel tot expressie komende genen werden gebruikt als input voor pathview. Een genvector met log2FC-waarden werd gekoppeld aan KEGG-pathways, met nadruk op hsa05323 (reumatoïde artritis pathway). Visualisaties werden gegenereerd waarin up- en downregulatie met kleur werd weergegeven.

pathview(
  gene.data = gene_vector,
  pathway.id = "hsa05323",
  species = "hsa",
  gene.idtype = "SYMBOL"
)

🧠 GO-analyse

Voor gene ontology (GO)-verrijkingsanalyse werd het pakket goseq gebruikt, inclusief correctie voor genlengtebias met behulp van genlengtes opgehaald via biomaRt. Enkel genen met significante expressieveranderingen (padj < 0.05 & |log2FC| > 1) werden meegenomen.

pwf <- nullp(all_genes, "hg38", "ensGene", bias.data = gene_lengths)
GO.wall <- goseq(pwf, "hg38", "ensGene")

📜 Analyseworkflow (stroomschema)

FASTQ (RA + Gezond)
     ↓
Mapping (Rsubread::align → GRCh38)
     ↓
BAM → Sorteren + indexeren
     ↓
featureCounts (GTF GRCh37.p13)
     ↓
DESeq2 → DE-analyse
     ↓
↓                    ↓
KEGG (pathview)     GO (goseq)

🧠 Documentatie

Het volledige script met commentaar en uitleg over elke stap is terug te vinden in scripts/RA_analysis_script.R. Alle parameters, bestandsnamen en filtercriteria zijn reproduceerbaar en consistent met bovenstaande beschrijving.



---


## 📊 Resultaten

🧬 Differentiële genexpressie

Na normalisatie en statistische analyse met DESeq2 werden in totaal ___ genen gevonden die significant differentieel tot expressie kwamen bij RA-patiënten ten opzichte van gezonde controlepersonen (padj < 0.05 en |log₂FC| > 1). Van deze genen waren er:

🔺 ___ genen up-gereguleerd

🔻 ___ genen down-gereguleerd

De genen met de sterkste expressieveranderingen waren onder andere [GENE1], [GENE2] en [GENE3].

📍 Volcano plot
In Figuur 1 is de volcano plot weergegeven. Op de x-as staat de log₂ fold change en op de y-as de -log₁₀ van de aangepaste p-waarde (padj). Significante genen zijn zichtbaar in kleur.

📁 Bestand: results/VolcanoPlot_RA_vs_Normal.png

Figuur 1. Volcano plot van differentiële genexpressie tussen RA en controle.


🔗 KEGG-pathwayanalyse
De genexpressieresultaten werden gekoppeld aan de KEGG-pathway "Rheumatoid Arthritis" (hsa05323). In Figuur 2 is te zien dat meerdere genen in deze pathway differentieel tot expressie komen, waaronder genen betrokken bij cytokinesignalering (TNF, IL1B, CXCL8).

📁 Bestand: results/pathview_hsa05323.png

Figuur 2. KEGG pathway visualisatie van hsa05323 met behulp van pathview. Rood: up-gereguleerd, groen: down-gereguleerd.

🧠 Gene Ontology (GO) analyse
De GO-enrichmentanalyse toonde significante oververtegenwoordiging van biologische processen gerelateerd aan het immuunsysteem. Enkele sterk verrijkte GO-termen zijn:

"immune response"

"leukocyte activation"

"response to cytokine"

📁 Bestand: results/GO_enrichment_plot.png

Figuur 3. Visualisatie van de top verrijkte GO-termen (biological process).

📋 Samenvatting

De analyse bevestigt dat RA geassocieerd is met grootschalige genexpressieveranderingen, met name in immuungerelateerde pathways en processen. Zowel de KEGG-pathway-analyse als GO-analyse ondersteunen de rol van ontsteking, cytokinesignalering en immuunactivatie in RA.


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




