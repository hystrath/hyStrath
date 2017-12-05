void createCellSet
(
    const polyMesh& mesh,
    const labelList& set,
    const word& setName
)
{
        // Load set to work
        autoPtr<topoSet> currentSetPtr(NULL);

        //- make a back up of the set
//         backup(mesh, setName, setName + "_old");

        currentSetPtr.reset
        (
            new cellSet
            (
                mesh,
                setName,
                set.size()
            )
        );

        topoSet& currentSet = currentSetPtr();

        // - filling set

        forAll(set, cell)
        {
            const label& cellI = set[cell];

            currentSet.insert(cellI);
        }

        currentSet.sync(mesh);

        Info<< "Writing " << currentSet.name()
            << " (size " << currentSet.size() << ") to "
            << currentSet.instance()/currentSet.local()
               /currentSet.name()
            << endl << endl;

        currentSet.write();
}
