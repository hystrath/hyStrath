# hyStrath

#### Hypersonic / Rarefied gas dynamics code developments released under license GPL-3.0 
#### The only platform to conjointly host open-source CFD and DSMC codes designed for atmospheric re-entry analysis

#### Includes:
+ *hyFoam* CFD solver (supersonic combusting flows)  
+ *hy2Foam* CFD solver (hypersonic reacting flows)  
+ *hy2MhdFoam* solver (hypersonic reacting flows with MHD flow control)  
+ **Coming soon: _ARC_**, a low computationally demanding 3-D Ablative Response Code  
+ *dsmcFoam+* code (direct simulation Monte Carlo)  
+ *pdFoam* code, a hybrid PIC-DSMC solver   

#### Please visit the [_hyStrath_ Wiki page](https://github.com/vincentcasseau/hyStrath/wiki)  

<br><br>

---  
## Compatibility, Maintenance, Installation and Sync

### Master branch  

#### Compatibility  
+ OF-v1706: [openfoam.com/releases/openfoam-v1706](https://www.openfoam.com/releases/openfoam-v1706)  

#### Installation  
```sh
git clone https://github.com/vincentcasseau/hyStrath.git --branch master --single-branch && cd hyStrath/  
./install-all.sh 8 > logInstall 2>&1 &
```  

where _8_ is the number of processors to be used during the installation.  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Compatibility,-Maintenance,-Installation-and-Sync)  

<br><br>

---  

## Publications

#### [How to cite our work](https://github.com/vincentcasseau/hyStrath/wiki/Publications#how-to-cite-our-work)  

#### Latest journal article:  
[*__dsmcFoam+__*] &nbsp; C. White _et al._, 03/2018: [dsmcFoam+: An OpenFOAM based direct simulation Monte Carlo solver](https://pure.strath.ac.uk/portal/files/81235392/White_etal_CPC_2017_an_OpenFOAM_based_direct_simulation_Monte_Carlo_solver.pdf)

#### Latest conference paper:  
[*__ARC__*] &nbsp; V. Renato _et al._, 09/2017: [Multi-dimensional Ablation and Thermal Response Program for Martian Entry
Analysis](https://strathprints.strath.ac.uk/62926)  

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/Publications)  


<br><br>

---  

## People & Contact

__GitHub coordinator:__ Dr Vincent Casseau  

__Contributors:__ Dr Vincent Casseau, Dr Daniel E.R. Espinoza, Dr Christopher J. Capon, Alexey Ryakhovskiy, Dr Jimmy-John O.E. Hoste, Dr Viola Renato, Dr Rodrigo C. Palharini, Dr Craig White, Dr Melrose Brown, Prof Russell R. Boyce, Dr Thomas J. Scanlon, Dr Richard E. Brown     

__External contributors:__ the Micro & Nano Flows Group, R.Tech Engineering   

#### [View more](https://github.com/vincentcasseau/hyStrath/wiki/People-and-Contact)  

#### [Contribute](https://github.com/vincentcasseau/hyStrath/wiki/Contributions)  


<br><br>

---  
### _hyStrath_ also features  
+ [**_blockMeshDG_**](https://openfoamwiki.net/index.php/Contrib_blockMeshDG) by Anton Kidess   
+ [**_makeAxialMesh_**](http://openfoamwiki.net/index.php/Contrib/MakeAxialMesh) by Bernhard Gschaider  
