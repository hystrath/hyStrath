# hyStrath

##### Hypersonic code developments in OpenFOAM released under license GPL-3.0
###### Hosting the *hyFoam* (supersonic combusting flows) and *hy2Foam* (hypersonic reacting flows) CFD solvers  
###### Hosting the *dsmcFoam+* (direct simulation Monte Carlo) solver  


---

**Please check out the [_hyStrath_ Wiki page](https://github.com/vincentcasseau/hyStrath/wiki)**   


---  
## Compatibility, Installation and Updates

###### Master branch (CFD and DSMC)  
+ OF-v1612+: https://www.openfoam.com/releases/openfoam-v1612+  

###### OF-2.4.0 branch (CFD only)   
+ OF-2.4.0-MNF: https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF  
+ OF-2.4.0: http://openfoam.org/download/2-4-0-ubuntu  
+ OF-2.3.0: http://openfoam.org/download/2-3-0-ubuntu  


<div class="paragraph"><p><br>
<br></p></div>


###### Installation  
OF-v1612+:  
+ git clone https://github.com/vincentcasseau/hyStrath.git --branch master --single-branch && cd hyStrath/   
+ CFD & DSMC: ./install-all.sh _#nCPUs_ > logInstall &
+ CFD: ./install-CFD.sh _#nCPUs_ > logInstall &
+ DSMC: ./install-DSMC.sh _#nCPUs_ > logInstall &   
 

OF-2.4.0-MNF, OF-2.4.0, OF-2.3.0:  
+ git clone https://github.com/vincentcasseau/hyStrath.git --branch OF-2.4.0 --single-branch && cd hyStrath/   
+ CFD: ./install.sh  

<div class="paragraph"><p><br>
<br></p></div>


###### Synchronisation
OF-v1612+:  
+ git pull origin master   
+ CFD & DSMC: ./resync-all.sh _#nCPUs_ > logSync &
+ CFD: ./resync-CFD.sh _#nCPUs_ > logSync &
+ DSMC: ./resync-DSMC.sh _#nCPUs_ > logSync & 


<div class="paragraph"><p><br>
<br></p></div>

---  
## Release history  
### 17 Feb 2018: 'Concordia' release, Master branch  
+ CFD upgrade  

#### 5 Dec 2017: 'Cairn' release, Master branch  
+ Introduction of the _dsmcFoam+_ solver  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Release-history)  


<div class="paragraph"><p><br>
<br></p></div>

---  

## Publications

#### Latest journal article:  
V. Casseau _et al._, 12/2016: [A Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows, Part Two: Multi-Dimensional Analysis](http://www.mdpi.com/2226-4310/3/4/45/html)  

#### Latest conference paper:  
J.-J. O.E. Hoste _et al._, 03/2017: [Numerical Modeling and Simulation of Supersonic Flows in Propulsion Systems by Open-Source Solvers](http://eprints.gla.ac.uk/140369/1/140369.pdf)  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Publications)  


<div class="paragraph"><p><br>
<br></p></div>

---  

## People and Contact

###### CFD  
Lead developer: Dr Vincent Casseau    
Contributors: Daniel E.R. Espinoza, Jimmy-John O.E. Hoste and Dr Thomas J. Scanlon              
Acknowledgements: Dr Rodrigo C. Palharini and Prof Richard E. Brown
  
   
###### DSMC        
Current developers: Dr Craig White and Dr Vincent Casseau    
Contributors: Daniel E.R. Espinoza and Dr Thomas J. Scanlon  
Acknowledgements: Dr Rodrigo C. Palharini  

**Enquiries: hy2Foam@gmail.com**  


<div class="paragraph"><p><br>
<br></p></div>

---  
### You may also be interested in:  
+ **_dsmcFoam+_** for OF-2.4.0-MNF   
the direct simulation Monte Carlo code from the Universities of Strathclyde and Glasgow  
hosted by Dr Craig White (https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF/tree/devel-craig)  

  
### _hyStrath_ also features:  
+ **_blockMeshDG_** by Akidess (https://openfoamwiki.net/index.php/Contrib_blockMeshDG)  
+ **_makeAxialMesh_** by Bernhard Gschaider (http://openfoamwiki.net/index.php/Contrib/MakeAxialMesh)


