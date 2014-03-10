#### This script is NOT intended to reproduce a particular analysis
#### Rather, it is the "sandbox" which was used to manually search
#### genes within gene sets. One by one, genes within
#### relevant gene sets were studied by using wikipedia, NCBI, and genecards.org
#### Based on these investigations, genes were further categorized
#### as is shown in Figure 3.
#### In the end, this script finishes with two variables: gtpasegenes and
#### metgenes. These varaibles are used in the Figure3.R script

load("./data/RData/cameraMod.RData")
load("./data/RData/cactus.RData")
load("./data/RData/cyberT.results.RData")
library(GO.db)
library(qvalue)
load("./data/annotation.RData/MsAnn.RData")
load("./data/annotation.RData/MsAnnGO.RData")

################################################################################
#                       Helper functions                                       #
################################################################################

# A versioin of the qvalue function that handles NA values by returning
# NA values
qvalue.na = function(pvals) {
  out = pvals*0
  qvals = qvalue(na.omit(pvals))$qvalue
  out[names(qvals)] = qvals
  return(out)
}

################################################################################
#                     Focusing on the GTPase categories                        #
################################################################################

### p and q values
ct = camera.result
qct = lapply(ct, function(x) apply(x,2,qvalue.na) )

### cycle through ontologies and samples to find enriched gtpase categories
i = 7 # change this
ont = "reactome" # and this to cycle through options
pv = sort(ct[[ont]][,i])
qv = qct[[ont]][names(pv),i]
sets = setnames = names(qv[ qv < 0.25 ])
if(ont %in% c("go","bp","cc","mf")) setnames = substr(Term(sets),1,50)
tab = data.frame( name=setnames, 
                  logpval=-log10(pv[sets]), 
                  qval=qv[sets] )
rownames(tab) = NULL
tab

### Manually save the enriched gtpase categories
gtpase.gos = c("GO:0005525", # GTP binding, A16 (q<0.25)
               "GO:0017048", # Rho GTPase binding, A16
               "GO:0003924", # GTPase activity, A16
               "GO:0005097", # Rab GTPase activator activity, A16
               "GO:0006184", # GTP catabaolic process, A16
               "GO:0048365", # Rac GTPase binding, A6-0.36 (qvalue)
               "GO:0030742", # GTP-dependent protein binding, A6-0.30
               "GO:0007265", # Ras protein signal transcduction, A6-0.28
               "GO:0032321", # positive regulation of Rho GTPase activity, A6-0.32
               "GO:0007264" ) # small GTPase signal transduction, A6-0.33

### get the Affymetrix probes in these ontology categories
mgis = unique( unlist( c( MsAnnGO$BP$MGI[gtpase.gos],
                          MsAnnGO$MF$MGI[gtpase.gos],
                          MsAnnGO$CC$MGI[gtpase.gos] ) ) )
affys = MsAnn$MGI$Affy[mgis]
affy.list = affys[ !sapply( affys, is.null ) ]
mgis = rep(names(affy.list),times=sapply(affy.list,length))
probes = unlist(affy.list)
symbols = as.character(sapply( MsAnn$MGI$Symbol[mgis], function(x) x[1] ))
genes = cbind(mgis,probes,symbols)
rownames(genes) = NULL

### Affy probe fold changes and p-values
library(gplots)
fc = cyberT.results$logFC[probes,]
rownames(fc) = symbols
pv = cyberT.results$pval[probes,]
qv = apply( pv, 2, function(x) p.adjust(x,"BH") )

# only look at genes with fold change>2 and q<0.01 (for sake of space)
fcmet = apply( abs(fc) > log2(1), 1, any ) # change these
qmet = apply( qv < 0.01, 1, any ) # to get a different list of genes

#### Now manually go through the genes included in the heatmap
#### See the code below the list of genes to generate this heatmap

################################################################################
#                            Manual list of gtpase genes                       #
################################################################################

### These "man.genes" are genes for which a categorization was not
### completely obvious
man.genes.2 = c("Plxnb1", # receptor for sema4d, helps activate rhoa
                "Cdc42ep2","Cdc42se1","Sema4d", # downstream of cdc42, actin
                "Grhl3", # transcriptonal RhoGEF19 regulator
                "Plk2", # phosphorylates to degrade rap inhibitor/activator and ras activator
                "Niacr1", # G protein coupled receptor
                "Ak4", # Regulates nucleotide composition
                "Inf2","Fmnl2","Ncf2","Diap3", # GTPase binding
                "Akap13", # adapter protein for Rho
                "Gmppb", # Enzyme GTP to GTP-mannose
                "Ehd4", # Dynamin gtpase domain, localizes to Rab5
                "Gch1", # GTP cyclohydrolase
                "Gem", # GTP binding protein
                "Pck1", # Enzyme which needs GTP to function
                "Dhcr24","Hmox1","Usp8","Jun","Insr","Erbb2","Srl", "Gtpbp8",# Unknown connection
                "Hbs1l","Gspt1", # GTP binding elongation factor
                "Gnl3", # Guanine nucleotide binding protein
                "Plce1", # Regulated by GTPases
                "Nme1", # GTP synthesis
                "Flna", # filamin, target of Rala
                "Guf1" # binds to mitochondrial ribosomes in gtp dependent manner
)
not.genes.2 = c("Ehd2") # binds ATP not GTP

ifn.gtpases = c("Gbp2","Gbp3","Iigp1","Tgtp1","Igtp","Irgm1","Irgm2")
small.gtpases = c("Rhov","Rhof","Rap2b","Arf2","Rab31","Rhob","Arl5b", "Rasl2-9",
                  "Ran","Rab27a", "Rab8b", "Rap1b", "Rras2", "Rac2", "Arl4a",
                  "Rasef", "Kras", "Rab10", "Rab17", "Rasl11b", "Rhod", "Rap2a",
                  "Rala", "Rnd1", "Rab40c")
gaps = c("Tbc1d24","Tbc1d9","Tbc1d2b","Tbc1d1","Iqgap2", "Grlf1", 
         "Tbc1d15","Tbc1d7")
gefs = c("Rapgef4","Rasgef1b","Rapgef3","Rapgef6","Vav3","Dock11","Dock5", 
         "Rapgef2","Net1", "Arhgef10l","Arhgef7", "Arhgef19","Dock8", "Arhgef2")
gproteins = c("Gna14", "Tgm2", "Gnal", "Gna11", "Gnaq","Gnas")
tubulin = c("Tuba4a","Tubb6","Tubb2b","Tuba1c","Tuba1b", "Tubb2a")
other.gtpases = c("Rem2","Gtpbp4","Sept9","Sept11")
gdis = c()

above2 = list( ifn.gtpases, small.gtpases, gaps, gefs, gproteins, tubulin, 
               other.gtpases, man.genes.2, not.genes.2, gdis)

#### Below FC of 2

man.genes = c("Eif2s3x","Rfxank","Usp6nl","Pkp4", "Hsp90aa1","Supt16h", "Lrrk1", #Unknown connection
              "Dvl1","Eif5b","Sdcbp","Vcl","Ulk1","Nin","Nkiras1","Klk1b4",
              "Hsp90ab1","Mapkapk5", "Mtl1","Nudt2", # unknown connection
              "Adss","Adssl1","Suclg1","Suclg2","Nkiras2","Pck2", # enzyme requiring GTP
              "Glud1", # inhibitid by GTP
              "Grb2", # recruits Rab5
              "Ehd1", # dynamin gtpase domain
              "Sar1a","Sar1b", # GTP dependent
              "Bcar3","Daam1","Rsu1","Pard6a","Diap1","Diap2","Pkn1","Rassf1", 
              "Eea1","Sgsm3", # GTPases signal transduction
              "Rtkn","G3bp2", # Scaffold for gtapse proteins
              "Ak3", # GTP:ATP phosphotransferase
              "Eif5","Mtif2", #GTPase activity
              "Thy1", # modulates Rho activity
              "Vps4a" # Rnd2 binding partner
)
not.genes = c() 

ifn.gtpases = c("Mx1")
small.gtpases = c("Nras","Cdc42","Rhoq","Rab24","Rab22a","Rhobtb3","Rab20","Rab43",
                  "Rab15","Rab11a","Rab7","Rnd3","Rhou","Rhoa","Rab4a","Rab19",
                  "Rabl3","Rab23","Rhot2","Arf4","Rab27b","Rhebl1","Rab14","Rap2c",
                  "Arf3","Arl5a","Arf6","Arl2","Arl6","Rac3","Arl15","Rhoc",
                  "Arl3","Rab7l1","Rras","Rab12","Rab2a","Rnd2","Rab33a","Rasd2",
                  "Arl5a","Rab44","Rab6a","Rab1","Arl8a","Rac1","Arl1",
                  "Rap1a","Rab33b","Rabl6","Arfrp1")
gaps = c("Rabgap1l","Tbc1d10a","Srgap2","Arhgap5","Arhgap8","Tbc1d4","Tbc1d12",
         "Tbc1d10b","Iqgap1","Rp2h","Tbc1d5","Tbc1d13","Tbc1d23","Tbc1d9b","Tbc1d17",
         "Agap3","Dlc1","Rabgap1","Tbc1d8b")
gefs = c("Dock6","Dock1","Dock7","Als2","Rasgrf2","Ralgps2","Rgl2",
         "Vav2","Rasgrp1","Rapgef5","Dock9","Ect2","Dock10","Rapgefl1","Sos1",
         "Sos2","Arhgef18","Rasgrp3","Arhgef9")
gdis = c("Rap1gds1","Ralgds","Gdi2")
gproteins = c("Gnaz","Gnai1","Gng2")
tubulin = c("Tubb5","Tubb4b","Tuba1a","Tuba8")
other.gtpases = c("Atl3","Gnl2","Sept7","Gpn1","Gfm1","Sept14","Trim23","Gimap1",
                  "Ift27","Gspt2","Eftud2","Dnm1l","Eftud1","Gnl1","Gimap8",
                  "Gimap6","Ola1","Dnm2","Rit1","Tube1","Sept4","Rem1",
                  "Gimap5","Gucy1a3","Atl1","Atl2","Mtg1","Spag1","Gfm2","Opa1")

below2 = list( ifn.gtpases, small.gtpases, gaps, gefs, gproteins, tubulin, 
               other.gtpases, man.genes, not.genes, gdis)

man.genes = c(unlist(below2),unlist(above2))

#### Now show genes left over 
fc2 = fc[ fcmet & qmet, ]
fc2 = fc2[ !(rownames(fc2) %in% c(man.genes,not.genes)), ] # helpful line to go through proteins one by one
#.2(fc2,col="bluered",symbreaks=TRUE,margins=c(5,5))

#### combine all the sig genes

names(below2) = c("ifn","small","gap","gef","G","tub","other","unknown","not","gdi")
names(above2) = c("ifn","small","gap","gef","G","tub","other","unknown","not","gdi")
sig.genes = mapply(c, below2, above2, SIMPLIFY=FALSE)
sig.genes$rho.binding = c("Inf2","Flna","Dock11","Akap13","Iqgap1",
                          "Diap1","Vcl","Diap3","Fmnl2","Dock11","Daam1",
                          "Rtkn","Vps4a","Dock9","Arhgef2","Pard6a",
                          "Fmnl2","Srgap1","Cit","Kif3b","Pfn2","Cdc42bpb",
                          "Lrrk2","Diap2","Dock10","Rtkn")
gtpasegenes = lapply(sig.genes,sort)
sapply( sig.genes, length )

################################################################################
#                  Focusing on the metabolic categories                        #
################################################################################

### p and q values
sc = cactus.result
qsc = sapply(sc, function(x) apply(x,2,qvalue.na) )

i = 6
ont = "bp"
pv = sort(sc[[ont]][,i])
qv = qsc[[ont]][names(pv),i]
sets = names(qv[ qv < 0.01 ])
if(ont %in% c("go","bp","cc","mf")) setnames = substr(Term(sets),1,50)
tab = data.frame( name=setnames, 
                  logpval=-log10(pv[sets]), 
                  qval=qv[sets] )
tab

gos = c("GO:0008203", # cholesterol metabolic process, A2-0.06
        "GO:0006695", # cholesterol biosynthetic process, A2-0.06
        "GO:0016126", # sterol biosynthetic process, A2-0.06
        "GO:0042953", # lipoprotein transport, A2-0.10
        "GO:0051384", # response to glucocorticoid stimulus, A2-0.11
        "GO:0008299", # isoprenoid biosynthetic process, A2-0.12
        "GO:0005179", # hormone activity, A2-0.26
        "GO:0050750", # low-density lipoprotein particle receptor binding, A2-0.26
        "GO:0005041", # low-density lipoprotein receptor activity, A2-0.26
        "GO:0034362", # low-density lipoprotein particle, A2-0.32
        "GO:0042627", # chylomicron, A2-0.32
        
        "GO:0051384", # response to glucocorticoid stimulus AB2-0.0003
        "GO:0008203", # cholesterol metabolic process, AB2-0.004
        "GO:0042953", # lipoprotein transport, AB2-0.05
        "GO:0008202", # steroid metabolic process, AB2-0.06
        "GO:0070328", # triglyceride homeostasis, AB2-0.06
        "GO:0006695", # cholesterol biosynthetic process, AB2-0.06
        "GO:0045723", # positive regulation of fatty acid biosynthetic pro, AB2-0.06
        "GO:0050995", # negative regulation of lipid catabolic process, AB2-0.09
        "GO:0008299", # isoprenoid biosynthetic process, AB2-0.10
        "GO:0043401", # steroid hormone mediated signaling pathway, AB2-0.14
        "GO:0046965", # retinoid X receptor binding, AB2-0.084117229
        "GO:0005504", # fatty acid binding, AB2-0.100820231
        "GO:0042974", # retinoic acid receptor binding, AB2-0.100820231
        "GO:0005041", # low-density lipoprotein receptor activity, AB2-0.121800141
        "GO:0035259", # glucocorticoid receptor binding, AB2-0.121800141
        "GO:0003707", # steroid hormone receptor activity, AB2-0.146974308
        "GO:0001972", # retinoic acid binding, AB2-0.161158112
        "GO:0070330", # aromatase activity, AB2-0.213738076
        "GO:0004364", # glutathione transferase activity, AB2-0.238179583
        "GO:0005184", # neuropeptide hormone activity, AB2-0.238179583
        "GO:0016614", # oxidoreductase activity, acting on CH-OH group of, AB2-0.238179583
        "GO:0004622", # lysophospholipase activity, AB2-0.256083763
        
        "GO:0004364", # glutathione transferase activity, A6-0.17
        "GO:0005496", # steroid binding, A6-0.17
        "GO:0016831", # carboxy-lyase activity, A6-0.17
        "GO:0004806", # triglyceride lipase activity, A6-0.19
        "GO:0070330", # aromatase activity, A6-0.20
        "GO:0004497", # monooxygenase activity, A6-0.20
        "GO:0015485", # cholesterol binding, A6-0.20
        "GO:0050253", # retinyl-palmitate esterase activity, A6-0.21
        "GO:0000062", # fatty-acyl-CoA binding, A6-0.23
        "GO:0034185", # apolipoprotein binding, A6-0.25
        "GO:0016817", # hydrolase activity, acting on acid anhydrides, A6-0.25
        "GO:0005041", # low-density lipoprotein receptor activity, A6-3.295734e-05
        "GO:0070330", # aromatase activity, A6-5.062006e-05
        "GO:0016810", # hydrolase activity, acting on carbon-nitrogen but, A6-3.097566e-04
        "GO:0034185", # apolipoprotein binding, A6-4.182055e-04
        "GO:0016705", # oxidoreductase activity, acting on paired donors, A6-5.458270e-04
        "GO:0030169", # low-density lipoprotein particle binding, A6-8.336624e-04
        "GO:0004497", # monooxygenase activity, A6-1.263385e-03
        "GO:0001972", # retinoic acid binding, A6-1.666315e-03
        "GO:0005496", # steroid binding, A6-1.691835e-03
        "GO:0004553", # hydrolase activity, hydrolyzing O-glycosyl compoun, A6-2.068417e-03
        "GO:0016788", # hydrolase activity, acting on ester bonds, A6-3.168354e-03
        "GO:0035259", # glucocorticoid receptor binding, A6-3.174952e-03
        "GO:0004029", # aldehyde dehydrogenase (NAD) activity, A6-3.296256e-03
        "GO:0004602", # glutathione peroxidase activity, A6-3.795437e-03
        "GO:0016818", # hydrolase activity, acting on acid anhydrides in, A6-5.445677e-03
        "GO:0004364", # glutathione transferase activity, A6-6.892754e-03        
        
        "GO:0005041", # low-density lipoprotein receptor activity, A6-0.30
        "GO:0070330", # aromatase activity, A6-0.31
        "GO:0016811", # hydrolase activity, acting on carbon-nitrogen (but, A6-0.31
        
        "GO:0019915", # lipid storage, A16-0.18
        #"GO:0006629") # lipid metabolic process, A16-0.28
        "GO:0006805", # xenobiotic metabolic process, A16-0.37
        "GO:0019395", # fatty acid oxidation, A16-0.38
        "GO:0006695", # cholesterol biosynthetic process, A16-0.38
        "GO:0006631", # fatty acid metabolic process, A16-0.38
        
        "GO:0004364", # glutathione transferase activity, A16-0.17
        "GO:0005496", # steroid binding, A16-0.17
        "GO:0016831", # carboxy-lyase activity, A16-0.17
        "GO:0004806", # triglyceride lipase activity, A16-0.19
        "GO:0070330", # aromatase activity, A16-0.20
        "GO:0004497", # monooxygenase activity, A16-0.20
        "GO:0015485", # cholesterol binding, A16-0.20
        "GO:0050253", # retinyl-palmitate esterase activity, A16-0.21
        "GO:0000062", # fatty-acyl-CoA binding, A16-0.23
        "GO:0034185", # apolipoprotein binding, A16-0.25
        "GO:0016817", # hydrolase activity, acting on acid anhydrides, A16-0.25
        "GO:0003707", # steroid hormone receptor activity, A16-0.25
        "GO:0031996", # thioesterase binding, A16-0.25
        "GO:0004602", # glutathione peroxidase activity, A16-0.27
        "GO:0004091") # carboxylesterase activity, A16-0.30

### get the Affymetrix probes in these ontology categories
gos = unique(gos)
mgis = unique( unlist( c( MsAnnGO$BP$MGI[gos],
                          MsAnnGO$MF$MGI[gos],
                          MsAnnGO$CC$MGI[gos] ) ) )
affys = MsAnn$MGI$Affy[mgis]
affy.list = affys[ !sapply( affys, is.null ) ]
mgis = rep(names(affy.list),times=sapply(affy.list,length))
probes = unlist(affy.list)
symbols = as.character(sapply( MsAnn$MGI$Symbol[mgis], function(x) x[1] ))
genes = cbind(mgis,probes,symbols)
rownames(genes) = NULL

### Affy probe fold changes and p-values
library(gplots)
fc = cyberT.results$logFC[probes,]
rownames(fc) = symbols
pv = cyberT.results$pval[probes,]
qv = apply( pv, 2, function(x) p.adjust(x,"BH") )

fcmet = apply( abs(fc) > log2(1), 1, any )
qmet = apply( qv < 0.01, 1, any )

################################################################################
#                            Manual list of metabolic genes                    #
################################################################################

man.genes2 = c( "Sorl1", # Binds LDL transports it into cells
                "Apoc2", # apolipoprotein
                "Apobec1", # edits Apob mRNA
                "Cubn", #Cotransporter role in lipoprotein, vitamin and iron metabolism
                "Zdhhc17", # maybe involved in endocytosis initiation
                "Lrpap1", # interacts with Lrp1
                "Mesdc2", # helps Lrp receptors
                "Hsp90b1","Dnaja1","Mmp13","Dkk1","Mesdc2", # LDL binding
                "Lrp8", # LDLR family
                "Cxcl16", # specifically binds to OxLDL
                "Olr1", # low density lipoprotein receptor
                "Lpl", # tryglyceride hydrolase + LDL uptake
                "O3far1", # omega fatty acid receptor
                
                "Dgat1","Dgat2", # triacylglyceral synthesis
                "Acox1","Acads","Acadm", # FA beta oxidation
                "Fabp2","Fabp6","Fabp4", # FA binding protein
                "Aadac","Lipg","Ces1d", # tryglyceride lipase activity 
                "Acer1","Acer3", # sphingolipid to FA
                "Plin2", # lipid storage
                "Hao2", # oxidation of hydroxyacids
                "Adipor2", # adiponectin receptor, lipid storage
                "Acsl4","Acsl5", # FA ligase
                "Acot11", # thioesterases, FA
                "Aasdh", # initial reaction in FA met
                "Ehhadh","Hadh","Acaa2","Ech1","Cpt1a","Hsd17b4","Echs1", # beta oxidation,
                "Slc7a2", # FA to FA ester
                "Acaa1b","Acaa1a", # beta oxidation pathway
                "Acsl3","Acsm3", # FA ligase
                "Lpin2","Acsf2", # FA met
                "Slc27a2","Slc27a4", # FA transporter
                "Lypla2", # hydrolyze FA
                "Fa2h", # FA hydroxylase, synthesis of FA containing sphingolipids
                "Aacs", # FA synthesis
                "Acot9","Acot10","Acot1","Acot2", # thioesterase
                "Mgll", # monoglyceride lipase
                
                "Stard4","Stard5", # Cholesterol binding
                "Cyp27a1", #  Catalyzes 1st step in oxidation of sterol intermediates
                "Lss", # conversion of (S)-2,3 oxidosqualene to lanosterol (cholest biosynthesis)
                "Hsd17b7" , # dehydrogenase, biosynthesis of cholesterol and steroids
                "Pmvk", # fifth reaction of the cholesterol biosynthetic pathway (mevalonate)
                "Fdps", # produces farnesyl pyrophosphate, intermediate in sterol biosynthesis (mevalonate)
                "Hmgcs1","Hmgcs2","Hmgcr", # (mevalonate)
                "Fdft1", #  first specific enzyme in cholesterol biosynthesis (mevalonate)
                "Idi1", # catalyzes the interconversion of isopentenyl diphosphate
                "Dhcr7", # 7-dehydrocholesterol to cholesterol (mevalonate)
                "Cyp51", # C14-demethylation of lanosterol
                "Prkaa1", # phosphorylates/stimulates catabolic enzymes
                "Sc4mol", # believed to function in cholesterol biosynthesis
                "Hpgd", # metabolism of prostaglandins
                "Gpam", # glycerolipid biosynthesis
                
                "Il1rn","Areg","Tnf","Il6","Fas", # glucocorticoid stimulus
                "Tat", # glucocorticoid stimulus
                "Cpn1", # protease, anaphylatoxin inactivator, inflam
                "Cdo1", #Initiates met paths relate pyruvate and sulfurate compounds
                "Bche", # esterase, detox, inactivate acetycholine
                "Pik3r1", # kinase, insulin interactions
                "Il1b", # neg reg of lipid catab
                "Niacr1", # nicotinic acid receptor
                "Hnf4g", # transcription factor
                "Cebp","Stat3", # glucocorticoid receptor binding
                "Xdh","Chdh","Dcxr","Gpd1l", # oxidorectase activity
                "Hdc","Amd2","Odc1","Pck1","Ddc","Uxs1","Pisd","Csad","Urod", # carboxy-lyase activity
                "Thbs1", # adhesive glycoprotein
                "Chi3l4","Chi3l3","Chi3l1","Hexb","Ganc","Idua","Hexb", # hydrolase activity
                "Snd1","Exo1","Pgap1","Aspa","Nt53", # hydrolase activity
                "Acvr1c", # activin receptor
                "Snca", # synuclein
                "Cryl1", # uronate pathway
                "Ghr", # GHR receptor
                "Pla2g15", # lipophospholipase
                "Gpx1", # glutathione peroxidase
                "Abhd3", # anhydrolase
                "Siae", # sialic acid acetylesterase
                
                "Adm","Gcg","Pyy","Prl2a1","Prl6a1","Hamp","Edn1","Edn2","Edn3",  # hormone activity
                "Ttr", # hormone/retinol transporter
                "Rdh1", # synthesis of retinaldehyde
                "Rora","Nr1h4","Rorc", # nuclear hormone receptor, retinoid related
                "Rarb","Ppara","Vdr", # RARbeta, PPARalpha, VDR
                "Nr6a1","Nr4a1","Nr5a2","Esrrg", # orphan nuclear receptor
                "Adra2a", # adranergic receptor
                "Smarca2", # VDR interaction
                
                "Iigp1","Irgm1", # gtpases
                "P4ha1","Ph4a2", # collagen chain folding
                "Dpysl3", # cytoskeleton
                "Egln3", # hydroxylates Hif1a/Hif2a
                "Vnn1","Upb1","Naaa","Nit2", # hydroxylases
                "Nudt7", # hydrolase
                
                "Lrat", # esterification of all-trans-retinol
                "Sult1a1", # sulfotransferase
                "Cyp2c55", # synthesis of 19-HETE, CAR PXR
                "Cyp2d22",
                "Ugt1a8","Adh5", # xenobiotics
                "Cyp4b1","Cyp2d9","Cyp2d26","Cyp2d10","Cyp2s1","Cyp3a25","Cyp2c40",
                "Ugt1a9","Ugt1a2","Ugt1a1","Ugt1a9","Ugt1a2","Ugt1a1","Cyp3a13", # aromatase
                "Gstm3","Mgst3","Gstk1","Mgst2","Gstt3","Gstm3","Gstm7","Gstm2",
                "Gstm1","Gstm6","Gstt1",
                "Ywhag","Tph1", # mono-oxygenase
                "Cyp20a1","Cyp2d34","Cyp2c68","Fmo1","Cyp4f14","Cyp4f15","Cyp4f13",
                "Cyp4v3","Fmo5","Cyp2w1", # mono-oxygenases
                "Ces1g","Ces1e", # carboxylesterases
                "Ces2b","Ces2a","Ces1f","Ces2c","Ces2e",
                "Aldh1a3","Aldh2","Aldh1b1", # detox aldehydes
                "Gpx2","Prdx6","Gpx4", # detox peroxides
                "Nceh1","Ugt1a6a", # xenobiotics
                "Por" # p450 enzymes
                
)

man.genes1 = c( "Apoa1","Apoe","Apob","Apoa2","Apoc3", # apoliprotein
                "Hdlbp", # HDL-binding protein
                "Vldlr","Ldlr","Lrp5","Lrp2","Lrp6", # LDL receptors
                "Abca1", # cholesteral efflux pump
                "Ldlrap1", # interacts with LDL receptor 
                "Lrp1", # endocytic receptor, lipid homeostasis
                "Pcsk9", # regulates cholesterol through ldlr
                "Lipc", # lipase, lipid uptake and hydrolysis
                "Lsr", # Binds chylomicrons, LDL and VLDL for uptake
                "Scarb1", # receptor many ligands, cholesterol, lipoproteins, etc
                
                "Nr1h3", # nuclear receptor interacts w/RXR
                "Srebf2", # TF, stimulates transcription of sterol-regulated genes
                "Scap", # Srebf chaperone
                "Cebpa","Mbtps1", # 1st step in the proteolytic activation of SREBP
                "Pparg","Rxra", # key upstream regulator
                "Asxl1", # RAR coactivator
                "Med1", # subunit of TR
                "Prmt2", # AR coactivator
                "Chd9", # Ppara coactivator
                
                "Nr4a3","Nr2c1","Nr2c2","Nr5a1","Nr2e3","Nr2c1", # Nuclear receptors
                "Nr3c2","Nr2f6","Nr3c2","Nr2f2","Nr1i3","Rxrb","Nr4a2","Ppard",
                "Ar","Nr3c1",
                "Nr1d2", # nuclear receptor
                "Ncoa1","Ncoa2", # coactivator steroid receptor
                "Nrip1","Ncor1", # interacts with nuclear receptor
                
                "Acox2","Acox3", # branched FA desaturation
                "Acadl","Acadvl", # FA metabolism
                "Hadhb","Hadha", # beta oxidation
                "Degs1", # FA desaturase
                "Acbd5", # binds medium-long chain acylCoA esters
                "Faah", # degrades bioactive FA
                "Gcdh", # glutaryl-CoA to crotonyl-CoA
                "Lypla1","Pnpla8","Pla2g4a","Lypla1","Pnpla6","Aspg","Enpp2", # lysophospholipase
                "Iah1", # lipase
                "Fasn", # fatty acid synthase
                "Lipa", # lipase
                "Mecr", # synthesis of FAs
                "Echdc2","Acot12","Agpat6", # FA met
                "Cpt1b", # transport of FAs
                "Lpin3", # regulates FA met
                "Cpt2", # LCFA met
                "Acsl1", # FA ligase
                "Crot","Crat", # beta oxidation
                "Acot4","Acot7", # FA thioesterase
                
                "Mttp", # tryglyceride transfer, lipoprotein assembly
                "Mvd","Mvk", # mevalonate decarboxylase and kinase
                "Soat2","Soat1", # Enzyme forms cholesterol esters
                "Sqle", # first oxygenation step in sterol biosynthesis
                "Cyp11a1", # enzyme cholesterol to pregnelone
                "Dhcr24", # Catalyzes the reduction sterol intermediates
                "Insig1", # control of cholesterol synthesis by controlling SCAP and HMGCR
                "Abcg1", # cholesterol and phospholipids transport
                "Nsdhl", # cholesterol biosynthesis
                "Prkaa2", # phosphorylates/stimulates catabolic enzymes
                "Cyb5r3", # desaturation and elongation of fattys, chol. biosynthesis, drug metab
                "Pdss1", # isoprenoid synthesis
                "Star", # role in hormone synthesis
                "Acbd3", # Steroid synth, golgi structure          
                
                
                "Mid1ip1", # lipid synthesis
                "Npc2", # cholesterol transport
                "Bscl2", # lipid storage
                
                
                "Pctp", # Catalyzes the transfer of phosphatidylcholine between membranes
                "Cebpa", # modulates leptin gene
                "Angptl3", # angiopoietin-like family
                "Cat", # converts ROS to water and oxygen
                "App", # Neural placticity/iron transfer
                "Cln8", # postulated function in lipid synthesis, transport, or sensing
                "Fech", # Catalyzes the ferrous insertion into protoporphyrin IX
                "Fdxr",  # initiates electron transport for cytochromes P450
                "Cftr",  # CF gene
                "Pappa", # protease cleaves IGF
                "Unc119b","Adam9","Ass1","Bcl2","Calcr","Agl","Map2k1", # glucocort stimulus
                "Pla2g7", # degreades Platelet activating factor
                "Slc37a4", # glucose homeostasis
                "Hnf4a","Mlxipl", # transcription factor
                "Gimap5", # gtpase
                "Syt1", # synaptic trafficking
                "Hmga1","Uimc1",
                "Cebpb","Ets2","Pias2", # glucocorticoid receptor binding
                "Ldhd","Agps", # oxidorectase activity
                "Paqr8","Paqr5","Serpina6","Pgrmc1", # steroid membrane receptor, glucort transport
                "Hsd11b2","Atp5o","Pgrmc2","Tor2a","Atp5o","Hsd17b10",
                "Ptch1", # Patched gene family
                "Hmgcl", # leucine degradation
                "Canx", # protein folding
                "Mapt", # microtubule associated protein
                "Plod3","Ogfod1","P4ha2","P4ha3","Amdhd2","Ino80","Dctd", #hydrolase activity
                "Smarca4","Dhx58","Brip1","Rad54l","Egln1","Ddx58","Supv3l1",
                "Ifih1","Plod1","Ydjc","Nadsyn1","Btd",
                "Ankra2", # facilitates endocytosis
                "Naga","Abhd10","Hexdc","Phka1","Chid1","Gbe1","Glb1","Phkb",
                "Man2c1","Man2b1","Hexa","Glb1", # hydrolase
                "Esd", # formaldehyde detoxification
                "Aoah","Supt6h","Olah","Vps29","Nt5e",
                "Pter", # phosphotriesterase
                "Acy3", # important role in deacetylating mercapturic acids
                "Chd4","Smarca5","Blm","Ddx5","Smarcal1", # hydrolase
                "Sirt2","Sirt3", # Sirtuins
                "Wdyhv1", # glutamine to glutamate
                "B4galnt1", # ganglioside biosynthesis
                "Mapk14", # a kinase
                "Cyb5", # hemoprotein
                "Nudt14", # UDP glucose
                "Arf6","Traf2","Traf4","Rac1","Traf3","Haus7","Traf6","Cdc42","Calm1", # Steroid hormone activity
                "Ppme1","Traf1","Pnpla7",
                
                "Aldh3a1", # detoxification acetaldehyde, met corticosteroids, biogenic amines, neurotransmitters, and lipid peroxidation.
                "Aldh1a1","Aldh1a7","Aldh3a2","Aldh7a1",
                
                "Ep300", # histone acetyltransferase, gluc stim, Hif1a activator
                
                "Hmgb1", # glucocorticoid stimulus
                "Plr6a1","Cck","Stc2","Retnla","Retnlb","Thpo","Lhb","Iapp",
                "Trh","Prl4a1","Sct","Prl7a1", # hormone activity
                "Nts","Vgf", # neuropeptide hormone
                "Umps","H47","Sgpl1","Pck2","Paics","Sgpl1","Ppcdc","Echdc1","Mlycd", # carboxylase activity
                
                "Dhrs4", # reduces all-trans-retinal and 9-cis-retinal
                "Cyp2c55", # synthesis of 19-HETE, CAR PXR
                "Cyp2d22","Cyp2b10","Cyp2c65","Cyp39a1", # steroid metabolism?
                "Sult1b1", # sulfotransferase
                "Srd5a1","Srd5a2", # testosterone to dihydrotestosterone
                "Igf2r","Cyp2j6","Cyp2j6", # aromatase activity
                "Gsta2","Gsta1","Gstm5","Gstz1","Gstm4","Gsto1","Gstt2","Hpgds", # Glutathione transferases
                "Ywhae","Ywhab","Cyp4a31","Cyp4a10","Cyp4a32","Fmo4","Cyp2d12", # monooxygenases
                "Cyp4a32","Cyp4a31","Cyp4a10","Cyp2c66","Cyp2j9","Cyp4f16", "Mical1",
                "Cyp2r1","Dohh","Mical2","Pam","Ces2g", # carboxylesterases
                "Akr1c12","Ugt2b1","Akr1c13", # xenobiotics
                
                "1110031I02Rik","AI464131","0610007P14Rik","Gm10639","Gm9705",
                "4932438A13Rik"
                
)

man.genes = c(man.genes2, man.genes1)

fc2 = fc[ fcmet & qmet, ]
fc2 = fc2[ !(rownames(fc2) %in% man.genes), ] # helpful line to go through proteins one by one
#heatmap.2(fc2,col="bluered",symbreaks=TRUE,margins=c(5,15))

# (1) uptake/export
uptake = c("Sorl1", "Apoc2", "Apobec1","Cubn","Zdhhc17","Lrpap1", 
           "Mesdc2","Hsp90b1","Dnaja1","Mmp13","Dkk1","Mesdc2",
           "Lrp8","Cxcl16","Olr1","Lpl","O3far1","Apoa1","Apoe",
           "Apob","Apoa2","Apoc3","Hdlbp","Vldlr","Ldlr","Lrp5",
           "Lrp2","Lrp6","Abca1","Ldlrap1","Lrp1","Pcsk9", 
           "Lipc","Lsr","Scarb1")

# (2) FA metabolism
FAmet = c("Dgat1","Dgat2","Acox1","Acads","Acadm","Fabp2","Fabp6",
          "Fabp4","Aadac","Lipg","Ces1d","Acer1","Acer3", 
          "Plin2","Hao2","Adipor2","Acsl4","Acsl5","Acot11", 
          "Aasdh","Ehhadh","Hadh","Acaa2","Ech1","Cpt1a","Hsd17b4",
          "Echs1","Slc7a2","Acaa1b","Acaa1a","Acsl3","Acsm3", 
          "Lpin2","Acsf2","Slc27a2","Slc27a4","Lypla2","Fa2h", 
          "Aacs","Acot9","Acot10","Acot1","Acot2","Mgll", 
          "Acox2","Acox3","Acadl","Acadvl","Hadhb","Hadha",
          "Degs1","Acbd5","Faah","Gcdh","Lypla1","Pnpla8","Pla2g4a",
          "Lypla1","Pnpla6","Aspg","Enpp2","Iah1","Fasn", 
          "Lipa","Mecr","Echdc2","Acot12","Agpat6","Cpt1b", 
          "Lpin3","Cpt2","Acsl1","Crot","Crat","Acot4","Acot7")

# (3) sterol/mevalonate
chol = c("Stard4","Stard5","Cyp27a1","Lss","Hsd17b7", 
         "Pmvk","Fdps","Hmgcs1","Hmgcs2","Hmgcr","Fdft1", 
         "Idi1","Dhcr7","Cyp51","Prkaa1","Sc4mol", 
         "Hpgd","Gpam","Mttp","Mvd","Mvk","Soat2","Soat1", 
         "Sqle","Cyp11a1","Dhcr24","Insig1","Abcg1","Nsdhl",
         "Prkaa2","Cyb5r3","Pdss1","Star","Acbd3",
         "Srebf2","Scap","Cebpa","Mbtps1")

# (4) detox enzymes
enzymes = c("Sult1a1","Cyp2c55","Cyp2d22","Ugt1a8","Adh5",
            "Cyp4b1","Cyp2d9","Cyp2d26","Cyp2d10","Cyp2s1","Cyp3a25","Cyp2c40",
            "Ugt1a9","Ugt1a2","Ugt1a1","Ugt1a9","Ugt1a2","Ugt1a1","Cyp3a13",
            "Gstm3","Mgst3","Gstk1","Mgst2","Gstt3","Gstm3","Gstm7","Gstm2",
            "Gstm1","Gstm6","Gstt1","Ywhag","Tph1","Cyp20a1","Cyp2d34",
            "Cyp2c68","Fmo1","Cyp4f14","Cyp4f15","Cyp4f13","Cyp4v3","Fmo5",
            "Cyp2w1","Ces1g","Ces1e","Ces2b","Ces2a","Ces1f","Ces2c","Ces2e",
            "Aldh1a3","Aldh2","Aldh1b1","Gpx2","Prdx6","Gpx4","Nceh1",
            "Ugt1a6a","Por","Dhrs4","Cyp2c55","Cyp2d22","Cyp2b10","Cyp2c65",
            "Cyp39a1","Sult1b1","Srd5a1","Srd5a2","Igf2r","Cyp2j6","Cyp2j6",
            "Gsta2","Gsta1","Gstm5","Gstz1","Gstm4","Gsto1","Gstt2","Hpgds",
            "Ywhae","Ywhab","Cyp4a31","Cyp4a10","Cyp4a32","Fmo4","Cyp2d12",
            "Cyp4a32","Cyp4a31","Cyp4a10","Cyp2c66","Cyp2j9","Cyp4f16", "Mical1",
            "Cyp2r1","Dohh","Mical2","Pam","Ces2g",
            "Akr1c12","Ugt2b1","Akr1c13")

# (5) nuclear receptors/regulators
NRs = c("Nr1h3","Srebf2","Scap","Cebpa","Mbtps1",
        "Pparg","Rxra","Asxl1","Med1","Prmt2","Chd9",
        "Nr4a3","Nr2c1","Nr2c2","Nr5a1","Nr2e3","Nr2c1",
        "Nr3c2","Nr2f6","Nr3c2","Nr2f2","Nr1i3","Rxrb","Nr4a2","Ppard",
        "Ar","Nr3c1","Nr1d2","Ncoa1","Ncoa2","Nrip1","Ncor1")

# (6) nocat
nocat = c("Pctp","Cebpa","Angptl3","Cat","App","Cln8", 
          "Fech","Fdxr","Cftr","Pappa","Unc119b","Adam9",
          "Ass1","Bcl2","Calcr","Agl","Map2k1","Pla2g7",
          "Slc37a4","Hnf4a","Mlxipl","Gimap5","Syt1",
          "Hmga1","Uimc1","Cebpb","Ets2","Pias2","Ldhd",
          "Agps","Paqr8","Paqr5","Serpina6","Pgrmc1",
          "Hsd11b2","Atp5o","Pgrmc2","Tor2a","Atp5o","Hsd17b10",
          "Ptch1","Hmgcl","Canx","Mapt","Plod3","Ogfod1",
          "P4ha2","P4ha3","Amdhd2","Ino80","Dctd","Smarca4",
          "Dhx58","Brip1","Rad54l","Egln1","Ddx58","Supv3l1",
          "Ifih1","Plod1","Ydjc","Nadsyn1","Btd","Ankra2",
          "Naga","Abhd10","Hexdc","Phka1","Chid1","Gbe1",
          "Glb1","Phkb","Man2c1","Man2b1","Hexa","Glb1",
          "Esd","Aoah","Supt6h","Olah","Vps29","Nt5e",
          "Pter","Acy3","Chd4","Smarca5","Blm","Ddx5","Smarcal1",
          "Sirt2","Sirt3","Wdyhv1","B4galnt1","Mapk14",
          "Cyb5","Nudt14","Arf6","Traf2","Traf4","Rac1","Traf3",
          "Haus7","Traf6","Cdc42","Calm1","Ppme1","Traf1","Pnpla7",
          "Aldh3a1","Aldh1a1","Aldh1a7","Aldh3a2","Aldh7a1",
          "Ep300","Hmgb1","Plr6a1","Cck","Stc2","Retnla",
          "Retnlb","Thpo","Lhb","Iapp","Trh","Prl4a1",
          "Sct","Prl7a1","Nts","Vgf","Umps","H47","Sgpl1",
          "Pck2","Paics","Sgpl1","Ppcdc","Echdc1","Mlycd",
          "Iigp1","Irgm1","P4ha1","Ph4a2","Dpysl3",
          "Egln3","Vnn1","Upb1","Naaa","Nit2","Nudt7")

metgenes = list(uptake=uptake,FAmet=FAmet,chol=chol,enzymes=enzymes,
                NRs=NRs,nocat=nocat)
metgenes = lapply(metgenes,sort)

################################################################################
#                       The end                                                #
################################################################################

rm( list=setdiff( ls(), c("metgenes","gtpasegenes") ) )






