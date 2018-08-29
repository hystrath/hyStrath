# hyStrath

#### Hypersonic / Rarefied gas dynamics code developments released under license GPL-3.0 
#### The only platform to conjointly host open-source CFD and DSMC codes designed for atmospheric re-entry analysis

#### Includes:
+ *hyFoam* CFD solver (supersonic combusting flows)  
+ *hy2Foam* CFD solver (hypersonic reacting flows)  
+ **Coming soon: _ARC_**, a low computationally demanding 3-D Ablative Response Code  
+ *dsmcFoam+* code (direct simulation Monte Carlo)  

#### Please visit the [_hyStrath_ Wiki page](https://github.com/vincentcasseau/hyStrath/wiki)  

<br><br>

---  
## Compatibility, Maintenance, Installation and Sync

### Master branch  

#### Compatibility  
+ OF-v1706: https://www.openfoam.com/releases/openfoam-v1706 

#### Installation  
```sh
git clone https://github.com/vincentcasseau/hyStrath.git --branch master --single-branch && cd hyStrath/  
./install-all.sh 8 > logInstall 2>&1 &
```  

where _8_ is the number of processors to be used during the installation.  

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

GitHub coordinator: Dr Vincent Casseau  

Contributors: Dr Vincent Casseau, Daniel E.R. Espinoza, Dr Jimmy-John O.E. Hoste, Dr Viola Renato, Dr Rodrigo C. Palharini, Dr Craig White, Dr Thomas J. Scanlon, Dr Richard E. Brown     
External contributors: [View more](https://github.com/vincentcasseau/hyStrath/wiki/Contributions)  


#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/People-and-Contact)  


<br><br>

---  
### _hyStrath_ also features  
+ **_blockMeshDG_** by Anton Kidess (https://openfoamwiki.net/index.php/Contrib_blockMeshDG)  
+ **_makeAxialMesh_** by Bernhard Gschaider (http://openfoamwiki.net/index.php/Contrib/MakeAxialMesh)
