# Hex sticker for PLMM
library(plmmr)
library(hexSticker)
library(showtext)
library(RColorBrewer)
library(corrplot)

RColorBrewer::brewer.pal(7, "Purples")
# "#F2F0F7" "#DADAEB" "#BCBDDC" "#9E9AC8" "#807DBA" "#6A51A3" "#4A1486"

png("patrick/corrplot_for_hex.png")
corrplot::corrplot(corr = cov2cor(pedigree$K[1:10, 1:10]),
                   col = RColorBrewer::brewer.pal(3, "Purples"),
                   # bg = "#EFEDF5",
                   tl.pos = "n",
                   cl.pos = "n")
dev.off()

font_add_google("Merriweather")
sticker("patrick/corrplot_for_hex.png",
        package="plmmr",
        p_color = "#3F007D",
        p_family = "Merriweather",
        p_size = 18,
        p_fontface = "bold",
        p_x = 0.83,
        p_y = 0.7,
        s_x=1,
        s_y=1,
        s_width = 0.6,
        s_height = 0.6,
        h_fill = "white",
        h_color = "#3F007D",
        h_size = 3,
        url = "pbreheny.github.io/plmmr/",
        u_color = "#54278F",
        u_size = 3,
        filename = "patrick/plmmr_hex_sticker.png")

