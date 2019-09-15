//---------------------------------------------------------
// Introduction of the MHD model
//---------------------------------------------------------
Info<< "\nReading MHD properties\n" << endl;
const dictionary mhdProperties =
(
    IFstream
    (
        fileName(thermo.lookup("mhdDictFile")).expand()
    )()
);

const word mhdType = mhdProperties.lookup("mhdModel");
const word conductivityType = mhdProperties.lookup("conductivityModel");
Info << "MHD model: "<< mhdType << endl;
Info << "Electrical conductivity model: " << conductivityType << endl;

Info<< "Creating MHD model\n" << endl;
autoPtr<mhd::mhdModel> mhd
(
    mhd::mhdModel::New
    (
        mhdProperties,
        thermo
    )
);
