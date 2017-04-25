# dropClust
## A package for automated quantification of multiplexed ddPCR data.

This R package aims at automatically quantifying the events of a multiplexed ddPCR reaction.
During a ddPCR run, each marker gene is fluorescently labeled with a combination of FAM and/or HEX fluorophore, giving it a unique footprint in the two-dimensional space represented by the intensities for each color channel. The position of each droplet within this space reveals, how many and, more importantly, which marker genes it contains. Thus, droplets that belong to the same marker cluster together. However, identifying and labelling these clusters correctly is not trivial, since one droplet can contain more than one marker. To accurately quantify the droplet count for each marker, it is crucial to both identify the clusters and label them correctly, based on their position. 
<p align="center">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387153/1b9928d2-29ca-11e7-97b3-ce67694eed0f.png"  
alt="Example B1"
width="400">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387163/20f98c90-29ca-11e7-8a9e-a7efaaf82fd9.png" 
alt="Example G1"
width="400">
</p>
