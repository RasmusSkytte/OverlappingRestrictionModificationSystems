# Combine panels into figures
# Figure 2
convert figures/Figure_2/fig2a.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" figures/Figure_2/fig2a.tif
convert figures/Figure_2/fig2b.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" figures/Figure_2/fig2b.tif
convert figures/Figure_2/fig2c.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "C" figures/Figure_2/fig2c.tif
#
convert figures/Figure_2/fig2a.tif figures/Figure_2/fig2b.tif figures/Figure_2/fig2c.tif +append figures/Figure_2/fig2.tif
#
rm figures/Figure_2/fig2a.tif
rm figures/Figure_2/fig2b.tif
rm figures/Figure_2/fig2c.tif
#
# Figure 3
convert Illustrations/fig3abc.png \( -append figures/Figure_3/fig3a.tif figures/Figure_3/fig3b.tif figures/Figure_3/fig3c.tif \) +append figures/Figure_3/fig3.tif
convert -flatten figures/Figure_3/fig3.tif figures/Figure_3/fig3.tif
#
convert figures/Figure_3/fig3.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +0+0    "A" figures/Figure_3/fig3.tif
convert figures/Figure_3/fig3.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +0+525  "B" figures/Figure_3/fig3.tif
convert figures/Figure_3/fig3.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +0+1050 "C" figures/Figure_3/fig3.tif
#
rm figures/Figure_3/fig3a.tif
rm figures/Figure_3/fig3b.tif
rm figures/Figure_3/fig3c.tif
#
# Figure 4
convert figures/Figure_4/fig4.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +30+0     "A" figures/Figure_4/fig4.tif
convert figures/Figure_4/fig4.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +670+0    "B" figures/Figure_4/fig4.tif
convert figures/Figure_4/fig4.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +1310+0   "C" figures/Figure_4/fig4.tif
convert figures/Figure_4/fig4.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +30+630   "D" figures/Figure_4/fig4.tif
convert figures/Figure_4/fig4.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +670+630  "E" figures/Figure_4/fig4.tif
convert figures/Figure_4/fig4.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +1310+630 "F" figures/Figure_4/fig4.tif
#
# Figure 5
cp figures/Figure_5/fig5_Bifidobacterium.tif figures/Figure_5/fig5a.tif # Aeromonas / Bifidobacterium / Clostridioides / Fusobacterium / Ralstonia / Yersinia
cp figures/Figure_5/fig5_Streptococcus.tif   figures/Figure_5/fig5b.tif # Bacillus / Bordetella / Enterococcus / Klebsiella / Pseudomonas / Staphylococcus / Streptococcus
cp figures/Figure_5/fig5_Paenibacillus.tif   figures/Figure_5/fig5c.tif # Chlamydia / Corynebacterium / Francisella / Paenibacillus / Salmonella?
#
convert figures/Figure_5/fig5a.tif figures/Figure_5/fig5b.tif figures/Figure_5/fig5c.tif +append figures/Figure_5/fig5abc.tif
convert figures/Figure_5/fig5d.tif figures/Figure_5/fig5e.tif figures/Figure_5/fig5f.tif +append figures/Figure_5/fig5def.tif
#
convert figures/Figure_5/fig5abc.tif figures/Figure_5/fig5def.tif -append figures/Figure_5/fig5.tif
#
rm figures/Figure_5/fig5abc.tif
rm figures/Figure_5/fig5def.tif
#
rm figures/Figure_5/fig5a.tif
rm figures/Figure_5/fig5b.tif
rm figures/Figure_5/fig5c.tif
rm figures/Figure_5/fig5d.tif
rm figures/Figure_5/fig5e.tif
rm figures/Figure_5/fig5f.tif
#
# Figure S1
convert figures/Figure_S1/figS1.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +50+40    "A" figures/Figure_S1/figS1.tif
convert figures/Figure_S1/figS1.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +1200+40  "B" figures/Figure_S1/figS1.tif
convert figures/Figure_S1/figS1.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +1200+520 "C" figures/Figure_S1/figS1.tif
#
# Figure S2
convert figures/Figure_S2/figS2.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +30+0   "A" figures/Figure_S2/figS2.tif
convert figures/Figure_S2/figS2.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +675+0  "B" figures/Figure_S2/figS2.tif
convert figures/Figure_S2/figS2.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +1320+0 "C" figures/Figure_S2/figS2.tif
#
# Figure S3
convert figures/Figure_S3/figS3a.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" figures/Figure_S3/figS3a.tif
convert figures/Figure_S3/figS3b.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" figures/Figure_S3/figS3b.tif
#
convert figures/Figure_S3/figS3a.tif figures/Figure_S3/figS3b.tif +append figures/Figure_S3/figS3.tif
#
rm figures/Figure_S3/figS3a.tif
rm figures/Figure_S3/figS3b.tif
#
# Figure S4
mkdir -p figures/Figure_S4 > /dev/null 2>&1
#
convert figures/Figure_5/fig5_Acinetobacter.tif      figures/Figure_5/fig5_Aeromonas.tif          figures/Figure_5/fig5_Bacillus.tif       +append figures/Figure_S4/panel1_row1.tif
convert figures/Figure_5/fig5_Bifidobacterium.tif    figures/Figure_5/fig5_Bordetella.tif         figures/Figure_5/fig5_Brucella.tif       +append figures/Figure_S4/panel1_row2.tif
convert figures/Figure_5/fig5_Burkholderia.tif       figures/Figure_5/fig5_Campylobacter.tif      figures/Figure_5/fig5_Chlamydia.tif      +append figures/Figure_S4/panel1_row3.tif
convert figures/Figure_5/fig5_Citrobacter.tif        figures/Figure_5/fig5_Clostridioides.tif     figures/Figure_5/fig5_Clostridium.tif    +append figures/Figure_S4/panel2_row1.tif
convert figures/Figure_5/fig5_Corynebacterium.tif    figures/Figure_5/fig5_Enterobacter.tif       figures/Figure_5/fig5_Enterococcus.tif   +append figures/Figure_S4/panel2_row2.tif
convert figures/Figure_5/fig5_Escherichia.tif        figures/Figure_5/fig5_Francisella.tif        figures/Figure_5/fig5_Fusobacterium.tif  +append figures/Figure_S4/panel2_row3.tif
convert figures/Figure_5/fig5_Haemophilus.tif        figures/Figure_5/fig5_Helicobacter.tif       figures/Figure_5/fig5_Klebsiella.tif     +append figures/Figure_S4/panel2_row4.tif
convert figures/Figure_5/fig5_Lactobacillus.tif      figures/Figure_5/fig5_Lactococcus.tif        figures/Figure_5/fig5_Legionella.tif     +append figures/Figure_S4/panel3_row1.tif
convert figures/Figure_5/fig5_Listeria.tif           figures/Figure_5/fig5_Mycobacterium.tif      figures/Figure_5/fig5_Mycoplasma.tif     +append figures/Figure_S4/panel3_row2.tif
convert figures/Figure_5/fig5_Neisseria.tif          figures/Figure_5/fig5_Paenibacillus.tif      figures/Figure_5/fig5_Pasteurella.tif    +append figures/Figure_S4/panel3_row3.tif
convert figures/Figure_5/fig5_Pseudomonas.tif        figures/Figure_5/fig5_Ralstonia.tif          figures/Figure_5/fig5_Rhizobium.tif      +append figures/Figure_S4/panel3_row4.tif
convert figures/Figure_5/fig5_Rickettsia.tif         figures/Figure_5/fig5_Salmonella.tif         figures/Figure_5/fig5_Serratia.tif       +append figures/Figure_S4/panel4_row1.tif
convert figures/Figure_5/fig5_Shewanella.tif         figures/Figure_5/fig5_Staphylococcus.tif     figures/Figure_5/fig5_Streptococcus.tif  +append figures/Figure_S4/panel4_row2.tif
convert figures/Figure_5/fig5_Streptomyces.tif       figures/Figure_5/fig5_Vibrio.tif             figures/Figure_5/fig5_Xanthomonas.tif    +append figures/Figure_S4/panel4_row3.tif
cp      figures/Figure_5/fig5_Yersinia.tif                                                                                                         figures/Figure_S4/panel4_row4.tif
#
convert figures/Figure_S4/panel1_row1.tif figures/Figure_S4/panel1_row2.tif figures/Figure_S4/panel1_row3.tif                                   -append figures/Figure_S4/figS4_1.tif
convert figures/Figure_S4/panel2_row1.tif figures/Figure_S4/panel2_row2.tif figures/Figure_S4/panel2_row3.tif figures/Figure_S4/panel2_row4.tif -append figures/Figure_S4/figS4_2.tif
convert figures/Figure_S4/panel3_row1.tif figures/Figure_S4/panel3_row2.tif figures/Figure_S4/panel3_row3.tif figures/Figure_S4/panel3_row4.tif -append figures/Figure_S4/figS4_3.tif
convert figures/Figure_S4/panel4_row1.tif figures/Figure_S4/panel4_row2.tif figures/Figure_S4/panel4_row3.tif figures/Figure_S4/panel4_row4.tif -append figures/Figure_S4/figS4_4.tif
#
rm figures/Figure_S4/panel1_row1.tif figures/Figure_S4/panel1_row2.tif figures/Figure_S4/panel1_row3.tif
rm figures/Figure_S4/panel2_row1.tif figures/Figure_S4/panel2_row2.tif figures/Figure_S4/panel2_row3.tif figures/Figure_S4/panel2_row4.tif
rm figures/Figure_S4/panel3_row1.tif figures/Figure_S4/panel3_row2.tif figures/Figure_S4/panel3_row3.tif figures/Figure_S4/panel3_row4.tif
rm figures/Figure_S4/panel4_row1.tif figures/Figure_S4/panel4_row2.tif figures/Figure_S4/panel4_row3.tif figures/Figure_S4/panel4_row4.tif
#
rm figures/Figure_5/fig5_Acinetobacter.tif      figures/Figure_5/fig5_Aeromonas.tif         figures/Figure_5/fig5_Bacillus.tif
rm figures/Figure_5/fig5_Bifidobacterium.tif    figures/Figure_5/fig5_Bordetella.tif        figures/Figure_5/fig5_Brucella.tif
rm figures/Figure_5/fig5_Burkholderia.tif       figures/Figure_5/fig5_Campylobacter.tif     figures/Figure_5/fig5_Chlamydia.tif
rm figures/Figure_5/fig5_Citrobacter.tif        figures/Figure_5/fig5_Clostridioides.tif    figures/Figure_5/fig5_Clostridium.tif
rm figures/Figure_5/fig5_Corynebacterium.tif    figures/Figure_5/fig5_Enterobacter.tif      figures/Figure_5/fig5_Enterococcus.tif
rm figures/Figure_5/fig5_Escherichia.tif        figures/Figure_5/fig5_Francisella.tif       figures/Figure_5/fig5_Fusobacterium.tif
rm figures/Figure_5/fig5_Haemophilus.tif        figures/Figure_5/fig5_Helicobacter.tif      figures/Figure_5/fig5_Klebsiella.tif
rm figures/Figure_5/fig5_Lactobacillus.tif      figures/Figure_5/fig5_Lactococcus.tif       figures/Figure_5/fig5_Legionella.tif
rm figures/Figure_5/fig5_Listeria.tif           figures/Figure_5/fig5_Mycobacterium.tif     figures/Figure_5/fig5_Mycoplasma.tif
rm figures/Figure_5/fig5_Neisseria.tif          figures/Figure_5/fig5_Paenibacillus.tif     figures/Figure_5/fig5_Pasteurella.tif
rm figures/Figure_5/fig5_Pseudomonas.tif        figures/Figure_5/fig5_Ralstonia.tif         figures/Figure_5/fig5_Rhizobium.tif
rm figures/Figure_5/fig5_Rickettsia.tif         figures/Figure_5/fig5_Salmonella.tif        figures/Figure_5/fig5_Serratia.tif
rm figures/Figure_5/fig5_Shewanella.tif         figures/Figure_5/fig5_Staphylococcus.tif    figures/Figure_5/fig5_Streptococcus.tif
rm figures/Figure_5/fig5_Streptomyces.tif       figures/Figure_5/fig5_Vibrio.tif            figures/Figure_5/fig5_Xanthomonas.tif
rm figures/Figure_5/fig5_Yersinia.tif
#
# Figure S5
convert figures/Figure_S5/figS5a.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" figures/Figure_S5/figS5a.tif
convert figures/Figure_S5/figS5b.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" figures/Figure_S5/figS5b.tif
#
convert figures/Figure_S5/figS5a.tif figures/Figure_S5/figS5b.tif +append figures/Figure_S5/figS5.tif
#
rm figures/Figure_S5/figS5a.tif
rm figures/Figure_S5/figS5b.tif
#
# Figure S6
convert figures/Figure_S6/figS6a.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" figures/Figure_S6/figS6a.tif
convert figures/Figure_S6/figS6b.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" figures/Figure_S6/figS6b.tif
convert figures/Figure_S6/figS6c.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "C" figures/Figure_S6/figS6c.tif
#
convert figures/Figure_S6/figS6a.tif figures/Figure_S6/figS6b.tif figures/Figure_S6/figS6c.tif +append figures/Figure_S6/figS6.tif
#
rm figures/Figure_S6/figS6a.tif
rm figures/Figure_S6/figS6b.tif
rm figures/Figure_S6/figS6c.tif
#
# Figure S7
convert figures/Figure_S7/figS7a.tif figures/Figure_S7/figS7b.tif +append figures/Figure_S7/figS7.tif
#
rm figures/Figure_S7/figS7a.tif
rm figures/Figure_S7/figS7b.tif
#
# Figure S8
convert figures/Figure_S8/figS8a.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "A" figures/Figure_S8/figS8a.tif
convert figures/Figure_S8/figS8b.tif   -gravity northwest -stroke none -fill black -pointsize 64 -annotate 0 "B" figures/Figure_S8/figS8b.tif
#
convert figures/Figure_S8/figS8a.tif figures/Figure_S8/figS8b.tif +append figures/Figure_S8/figS8.tif
#
rm figures/Figure_S8/figS8a.tif
rm figures/Figure_S8/figS8b.tif
#
# Figure S9
convert figures/Figure_S9/figS9_Escherichia_10_percent.tif   figures/Figure_S9/figS9_Escherichia_20_percent.tif     figures/Figure_S9/figS9_Escherichia_30_percent.tif     +append figures/Figure_S9/row1.tif
convert figures/Figure_S9/figS9_Escherichia_40_percent.tif   figures/Figure_S9/figS9_Escherichia_50_percent.tif     figures/Figure_S9/figS9_Escherichia_60_percent.tif     +append figures/Figure_S9/row2.tif
convert figures/Figure_S9/figS9_Escherichia_70_percent.tif   figures/Figure_S9/figS9_Escherichia_80_percent.tif     figures/Figure_S9/figS9_Escherichia_90_percent.tif     +append figures/Figure_S9/row3.tif
mv      figures/Figure_S9/figS9_Escherichia_100_percent.tif  figures/Figure_S9/row4.tif
#
convert figures/Figure_S9/row1.tif figures/Figure_S9/row2.tif figures/Figure_S9/row3.tif figures/Figure_S9/row4.tif -append figures/Figure_S9/figS9.tif
#
rm figures/Figure_S9/row1.tif figures/Figure_S9/row2.tif figures/Figure_S9/row3.tif figures/Figure_S9/row4.tif
#
rm figures/Figure_S9/figS9_Escherichia_10_percent.tif   figures/Figure_S9/figS9_Escherichia_20_percent.tif     figures/Figure_S9/figS9_Escherichia_30_percent.tif
rm figures/Figure_S9/figS9_Escherichia_40_percent.tif   figures/Figure_S9/figS9_Escherichia_50_percent.tif     figures/Figure_S9/figS9_Escherichia_60_percent.tif
rm figures/Figure_S9/figS9_Escherichia_70_percent.tif   figures/Figure_S9/figS9_Escherichia_80_percent.tif     figures/Figure_S9/figS9_Escherichia_90_percent.tif
#
# Figure S10
convert figures/Figure_S10/figS10bc.tif Illustrations/figS10b_insert.png -geometry +485+40  -composite figures/Figure_S10/figS10bc.tif
convert figures/Figure_S10/figS10bc.tif Illustrations/figS10c_insert.png -geometry +140+350 -composite figures/Figure_S10/figS10bc.tif
#
convert Illustrations/figS10a.png figures/Figure_S10/figS10bc.tif +append figures/Figure_S10/figS10.tif
#
convert -flatten figures/Figure_S10/figS10.tif figures/Figure_S10/figS10.tif
#
convert figures/Figure_S10/figS10.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +0+0     "A" figures/Figure_S10/figS10.tif
convert figures/Figure_S10/figS10.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +630+0   "B" figures/Figure_S10/figS10.tif
convert figures/Figure_S10/figS10.tif -gravity northwest -stroke none -fill black -pointsize 64 -annotate +630+310 "C" figures/Figure_S10/figS10.tif
#
rm figures/Figure_S10/figS10bc.tif
#
# Make png copy of tif files
find figures -name '*.tif' -print0 | xargs -0 -r mogrify -format png