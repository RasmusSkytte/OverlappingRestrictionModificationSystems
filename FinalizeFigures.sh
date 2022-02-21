# Overlay the illustrations
convert figures/Figure_5/fig5.png Illustrations/fig5.png -geometry +100+40 -composite figures/Figure_5/fig5.png

convert figures/Figure_S9/figs9.png Illustrations/figS9a.png -composite figures/Figure_S9/figs9.png
convert figures/Figure_S9/figs9.png Illustrations/figS9b_insert.png -geometry +1950+100 -composite figures/Figure_S9/figs9.png
convert figures/Figure_S9/figs9.png Illustrations/figS9c_insert.png -geometry +1600+850 -composite figures/Figure_S9/figs9.png