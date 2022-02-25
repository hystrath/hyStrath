//---------------------------------------------------------
// Introduction of the MHD model
//---------------------------------------------------------

dictionary mhdProperties = dictionary::null;
autoPtr<mhd::mhdModel> mhd(NULL);

if (thermo.found("mhdDictFile"))
{
    if (isFile(fileName(thermo.lookup("mhdDictFile")).expand()))
    {
        Info<< "Reading MHD properties" << endl;
        mhdProperties =
        (
            IFstream
            (
                fileName(thermo.lookup("mhdDictFile")).expand()
            )()
        );
        
        Info<< "Creating MHD model" << endl;
        mhd = mhd::mhdModel::New(mhdProperties, thermo);
    }
}

if (mhd.empty())
{
    Info<< "Creating MHD model" << endl;
    mhd = mhd::mhdModel::New(thermo);
}

mhd->update(U);
