gtm
===

Generative topographic map, reboot of Ondrej Such's R package

GTM is a cool technique; I like it a little better than SOM.
It's a shame there are none in CRAN, so I figure I'll bring this
one up to date and resubmit it.

Updates of original package: 
1) renamed the functions to get rid of _'s
2) commented out an apparently extraneneous gc()
3) added (necessary for R-3.x) NAMESPACE file
4) rewrote distance measure for speed (using proxy package)

ToDo:
0) Remove extraneous for loops in gtm.resp6, gtm.trn, add test scripts
1) Add some visualizations/OO classes; potential help in:
    https://github.com/geoss/som_visualization_r	
    http://stackoverflow.com/questions/19858729/r-package-kohonen-how-to-plot-hexagons-instead-of-circles-as-in-matlab-som-too
    http://nbremer.blogspot.nl/2013/07/on-creation-of-extended-self-organizing.html
    http://nbremer.blogspot.nl/2013/11/how-to-create-hexagonal-heatmap-in-r.html
    Also check out hexbin: http://cran.r-project.org/web/packages/hexbin/index.html
2) Possible heirarchical representations, aka
   http://eprints.aston.ac.uk/1314/1/IEEE_Transactions_24_%285%29.pdf
3) Temporal extension? 
   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.129.4734
