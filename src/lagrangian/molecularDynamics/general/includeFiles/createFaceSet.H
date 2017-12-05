void createFaceSet
(
    const polyMesh& mesh,
    const labelList& set,
    const word& setName
)
{
//     if(set.size() > 0)
//     {
        // Load set to work
        autoPtr<topoSet> currentSetPtr(NULL);

        //- make a back up of the set
//         backup(mesh, setName, setName+"_old");

        currentSetPtr.reset
        (
            new faceSet
            (
                mesh,
                setName,
                set.size()
            )
        );

        topoSet& currentSet = currentSetPtr();

        // - filling set
        forAll(set, face)
        {
            const label& faceI = set[face];

            currentSet.insert(faceI);
        }

        currentSet.sync(mesh); // removed this cause it was messing about with parallisation

        Info<< "Writing " << currentSet.name()
            << " (size " << currentSet.size() << ") to "
            << currentSet.instance()/currentSet.local()
               /currentSet.name()
            << endl << endl;

        currentSet.write();
//     }
}


// try to read in sets by the following code 

//         autoPtr<topoSet> currentSetPtr(NULL);
//         IOobject::readOption r;
// 
//         r = IOobject::MUST_READ;
// 
//         currentSetPtr.reset
//         (
//             new faceSet
//             (
//                 mesh,
//                 setName,
//                 r
//             )
//         );
