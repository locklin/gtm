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
1) Add some visualizations/OO classes
2) Possible heirarchical representations, aka
   http://eprints.aston.ac.uk/1314/1/IEEE_Transactions_24_%285%29.pdf
3) Temporal extension? 
   http://citeseerx.ist.psu.edu/viewdoc/summary?doi=10.1.1.129.4734
