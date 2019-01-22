//---------------------------------------------------------
// Introduction of the MHD model
//---------------------------------------------------------
const dictionary mhdProperties =
(
    IFstream
    (
        fileName(thermo.lookup("mhdDictFile")).expand()
    )()
);
const word mhdType = mhdProperties.lookup("mhdModel");
const word conductivityType = mhdProperties.lookup("conductivityModel");
Info << "Mhd model: "<< mhdType << endl;
Info << "Conductivity model: " << conductivityType << endl;

Info<< "Creating MHD model\n" << endl;
autoPtr<mhd::mhdModel> mhd
(
    mhd::mhdModel::New
    (
        mhdProperties,
        thermo
    )
);
