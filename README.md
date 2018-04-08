# hyStrath

#### Hypersonics / Rarefied gas dynamics code developments released under license GPL-3.0 
#### The only platform to conjointly host open-source CFD and DSMC codes designed for atmospheric re-entry analysis

#### Includes:
+ *hyFoam* CFD solver (supersonic combusting flows)  
+ *hy2Foam* CFD solver (hypersonic reacting flows)  
+ *dsmcFoam+* code (direct simulation Monte Carlo)  

#### Please visit the [_hyStrath_ Wiki page](https://github.com/vincentcasseau/hyStrath/wiki)  

<br><br>

---  
## Compatibility, Maintenance, Installation and Sync

### Master branch (CFD and DSMC)  

#### Compatibility  
+ OF-v1612+: https://www.openfoam.com/releases/openfoam-v1612+ 

#### Installation  
```sh
git clone https://github.com/vincentcasseau/hyStrath.git --branch master --single-branch && cd hyStrath/
```   
+ CFD & DSMC:
```sh 
./install-all.sh nCPUs > logInstall &
```

where _nCPUs_ is the number of processors to be used during the installation.  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Compatibility,-Maintenance,-Installation-and-Sync)  

<br><br>

---  
## Release history  
#### 17 Feb 2018: 'Concordia' release, Master branch  
+ CFD upgrade  
+ DSMC update 

##### 5 Dec 2017: 'Cairn' release, Master branch  
+ Introduction of the _dsmcFoam+_ solver  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Release-history)  


<br><br>

---  

## Publications

#### Latest journal article:  
V. Casseau _et al._, 12/2016: [A Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows, Part Two: Multi-Dimensional Analysis](http://www.mdpi.com/2226-4310/3/4/45/html)  

#### Latest conference paper:  
J.-J. O.E. Hoste _et al._, 03/2017: [Numerical Modeling and Simulation of Supersonic Flows in Propulsion Systems by Open-Source Solvers](http://eprints.gla.ac.uk/140369/1/140369.pdf)  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Publications)  


<br><br>

---  

## People and Contact

#### _hyStrath_ platform  
GitHub coordinator: Dr Vincent Casseau  

#### CFD  
Lead developer: Dr Vincent Casseau    
Contributors: Daniel E.R. Espinoza, Jimmy-John O.E. Hoste and Dr Thomas J. Scanlon              
Acknowledgements: Dr Rodrigo C. Palharini and Prof Richard E. Brown
     
#### DSMC        
Current developers: Dr Craig White and Dr Vincent Casseau    
Contributors: Daniel E.R. Espinoza and Dr Thomas J. Scanlon  
Acknowledgements: Dr Rodrigo C. Palharini  

**Enquiries: hy2Foam@gmail.com**  
**Consultancy services: [MTS-CFD](https://www.mts-cfd.com/)**  


<br><br>

---  
### _hyStrath_ also features  
+ **_dsmcFoam+_** for OF-2.4.0-MNF   
the direct simulation Monte Carlo code from the Universities of Strathclyde and Glasgow  
hosted by Dr Craig White (https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF/tree/devel-craig)  
+ **_blockMeshDG_** by Akidess (https://openfoamwiki.net/index.php/Contrib_blockMeshDG)  
+ **_makeAxialMesh_** by Bernhard Gschaider (http://openfoamwiki.net/index.php/Contrib/MakeAxialMesh)
