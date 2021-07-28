# Chemical family
tmp = c("Organochlorines", "Organophosphate metabolites", "Pyrethroids",
        "Pyrethroid metabolites", "PCBs", "PBDEs", "Acidic herbicides",
        "Anilinopyrimidines", "Azoles", "Benzamides", "Carbamates", "Carboxamides",
        "Neonicotinoids", "Oxadiazines", "Phenylpyrazoles", "Strobilurins",
        "Triazines/Triazinones/Diazines", "Subtituted ureas", "Dinitroanilines",
        "Other")
annot.colours = colorspace::lighten(colorRampPalette(RColorBrewer::brewer.pal(n=12,name='Paired'))(length(tmp)), amount=0.2)
names(annot.colours) = tmp

# Batches
filepaths = c("Luxembourg", "France", "GrandeSynthe", "Pooled3", "Pooled2")
batches = c("Luxembourg", "France", "Grande-Synthe", "Pooled", "LUX and GS")
batch.colours = c("tomato", "royalblue", "forestgreen", "mediumpurple", "darkgoldenrod")
names(batch.colours) = c("LUX","FRA","GS","Pooled3","Pooled2")

# Family
family.colours = c(colorRampPalette(RColorBrewer::brewer.pal(n=11,name='Spectral'))(66), "grey")
names(family.colours) = c(paste0("F",1:35), paste0("G",1:13), paste0("L",1:18), "Isolated")

# Region
area = read.csv("../Dictionaries/French_area_codes.csv")
region.colours = colorRampPalette(RColorBrewer::brewer.pal(n=9,name='Set1'))(length(unique(area$Region))+1)
names(region.colours) = c(unique(area$Region),"Luxembourg")

# Department
depart.colours = colorRampPalette(RColorBrewer::brewer.pal(n=9,name='Set1'))(length(unique(area$Department))+1)
names(depart.colours) = c(unique(area$Department),"Luxembourg")

# # Check colours
# image(1:length(depart.colours), 1, as.matrix(1:length(depart.colours)),
#       col=depart.colours,
#       xlab="", ylab = "", xaxt = "n", yaxt = "n", bty = "n")

