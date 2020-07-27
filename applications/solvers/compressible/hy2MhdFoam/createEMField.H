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

Info<< "Creating MHD model\n" << endl;
autoPtr<mhd::mhdModel> mhd
(
    mhd::mhdModel::New
    (
        mhdProperties,
        thermo
    )
);
