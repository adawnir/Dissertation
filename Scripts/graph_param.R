# Chemical family colours
chem_family = c("Organochlorines", "Organophosphate metabolites", "Pyrethroids",
                "Pyrethroid metabolites", "PCBs", "PBDEs", "Acidic herbicides",
                "Anilinopyrimidines", "Azoles", "Benzamides", "Carbamates", "Carboxamides",
                "Neonicotinoids", "Oxadiazines", "Phenylpyrazoles", "Strobilurins",
                "Triazines/Triazinones/Diazines", "Subtituted ureas", "Dinitroanilines",
                "Other")
annot.colours = lighten(colorRampPalette(RColorBrewer::brewer.pal(n=12,name='Paired'))(length(chem_family)), amount=0.2)
names(annot.colours) = chem_family

filepaths = c("Luxembourg", "France", "GrandeSynthe", "Pooled3", "Pooled2")
batches = c("Luxembourg", "France", "Grande-Synthe", "Pooled (LUX/FRA/GS)", "Pooled (LUX/GS)")
batch.colours = c("tomato", "royalblue", "forestgreen", "navy", "darkgoldenrod")