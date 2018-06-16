# hyStrath

#### Hypersonic / Rarefied gas dynamics code developments released under license GPL-3.0 
#### The only platform to conjointly host open-source CFD and DSMC codes designed for atmospheric re-entry analysis

#### Includes:
+ *hyFoam* CFD solver (supersonic combusting flows)  
+ *hy2Foam* CFD solver (hypersonic reacting flows)  
+ **Coming soon: ARC**, a low computationally demanding 3-D Ablative Response Code  
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
./install-all.sh nCPUs > logInstall 2>&1 &
```

where _nCPUs_ is the number of processors to be used during the installation.  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Compatibility,-Maintenance,-Installation-and-Sync)  

<br><br>

---  
## Release history  
#### 17 Feb 2018: 'Concordia' release, Master branch  
+ CFD upgrade to OpenFOAM v1612+   
+ DSMC update (the new features are listed [here](https://github.com/vincentcasseau/hyStrath/wiki)) 

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Release-history)  


<br><br>

---  

## Publications

#### Latest journal article:  
[*__hy2Foam__*] V. Casseau _et al._, 12/2016: [A Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows, Part Two: Multi-Dimensional Analysis](http://www.mdpi.com/2226-4310/3/4/45/html)  

#### Latest conference paper:  
[*__ARC__*] V. Renato _et al._, 09/2017: [Multi-dimensional Ablation and Thermal Response Program for Martian Entry
Analysis](https://strathprints.strath.ac.uk/62926)  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Publications)  


<br><br>

---  

## People & Contact

GitHub coordinator: Dr Vincent Casseau <a style="text-decoration: none" href="https://uk.linkedin.com/in/vincentcasseau" target="_blank"><img src="https://i2.wp.com/poxse.com/wp-content/uploads/2016/01/linkedin-logo.jpg?ssl=1" alt="LinkedIn profile" width="15"></a>
<a style="text-decoration: none" href="https://www.researchgate.net/profile/Vincent_Casseau" target="_blank"><img src="https://www.wur.nl/upload_mm/2/8/5/9f59698c-e156-4f33-9520-405cb7f4d9c6_researchgate_56f72ad6_490x330.png" alt="ResearchGate profile" width="22"></a>   

Lead CFD developer: Dr Vincent Casseau    
Current DSMC developer: Dr Vincent Casseau    
Lead ARC developer: Dr Viola Renato  

Contributors: Daniel E.R. Espinoza, Jimmy-John O.E. Hoste, Dr Craig White and Dr Thomas J. Scanlon    
External contributors: [View more](https://github.com/vincentcasseau/hyStrath/wiki/Contributions)  

**Enquiries: hystrath@gmail.com**   

#### [View all](https://github.com/vincentcasseau/hyStrath/wiki/People-and-Contact)  


<br><br>

---  
### _hyStrath_ also features  
+ **_dsmcFoam+_** for OF-2.4.0-MNF   
the direct simulation Monte Carlo code from the Universities of Strathclyde and Glasgow  
hosted by Dr Craig White (https://github.com/MicroNanoFlows/OpenFOAM-2.4.0-MNF/tree/devel-craig)  
+ **_blockMeshDG_** by Akidess (https://openfoamwiki.net/index.php/Contrib_blockMeshDG)  
+ **_makeAxialMesh_** by Bernhard Gschaider (http://openfoamwiki.net/index.php/Contrib/MakeAxialMesh)
