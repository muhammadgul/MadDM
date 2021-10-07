# MadDM
## install madgraph first but check the latest version first
 [Here](https://github.com/IPNL-CMS/HTTMadgraphDocumentation) are some important tips.

###### Install madgraph from browser, 
[Download](http://madgraph.phys.ucl.ac.be/) from browser.
###### Or from command line
```
wget http://madgraph.physics.illinois.edu/Downloads/MG5_aMC_v2.6.4.tar.gz
tar xf MG5_aMC_v2.6.4.tar.gz
rm MG5_aMC_v2.6.4.tar.gz
```
## MadDM installation
###### Using [paper](https://arxiv.org/pdf/1804.00044.pdf) to install maddm
###### Install the related [software](http://johannesbuchner.github.io/PyMultiNest/install.html) before installation of MadDM
###### Include LD_LIBRARY_PATH in **.bashrc** so that it will automatically included when open a new terminal
```export LD_LIBRARY_PATH=$PWD/MultiNest/lib:$PWD/cuba/directory/:$LD_LIBRARY_PATH```
./bin/mg5_aMC
install maddm

# run maddm
./bin/maddm.py
# example run:
import model DMsimp_s_spin0_MD,
define darkmatter ~xd,
generate relic_density
add direct_detection
add indirect_detection
output Example_1
launch Example_1
set sigmav_method inclusive
indirect=flux_source
set indirect_flux_source_method PPPC4DMID
set Mxd 1000
set My0 500

