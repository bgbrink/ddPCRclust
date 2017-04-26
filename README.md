# dropClust
## A package for automated quantification of multiplexed ddPCR data.

This R package was designed to automatically quantifying the events of a multiplexed ddPCR reaction.
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

For robustness, dropClust incorporates three independent clustering algorithms: the flowDensity algorithm, published in 2012 by M. Jafar Taghiyar and Mehrnoush Malek, the clustering algorithm SamSPECTRAL, a version of spectral clustering adapted to flow cytometry data and the flowPeaks package, developed by Yongchao Ge and Stuart C. Sealfon. The results are combined into a cluster ensemble. The agreement between the three approaches provides a measure of confidence for clustering results.
<p align="center">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387160/1e5eea02-29ca-11e7-871b-e2e3cd2639ec.png"  
alt="Result B1"
width="400">
<img 
src="https://cloud.githubusercontent.com/assets/11661112/25387164/224c4876-29ca-11e7-9f8f-557b7e515f0a.png" 
alt="Example G1"
width="400">
</p>

## Usage
Lorem ipsum

## License
    dropClust
    Copyright (C) 2017  Benedikt G. Brink

    This program is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License
    as published by the Free Software Foundation; either version 2
    of the License, or (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
