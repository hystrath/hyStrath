/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2016-2021 hyStrath
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of hyStrath, a derivative work of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "Linear.H"
#include "pointMesh.H"
#include "pointField.H"
#include "oppositeFace.H"
#include "coupledPointPatchField.H"

// * * * * * * * * * * * * * Static Member Functions * * * * * * * * * * * * //

template<class Type>
Foam::autoPtr<labelList> Foam::AveragingMethods::Linear<Type>::size
(
    const fvMesh& mesh
)
{
    autoPtr<labelList> s(new labelList(2));
    s()[0] = mesh.nCells();
    s()[1] = mesh.nPoints();
    return s;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Linear<Type>::Linear
(
    const IOobject& io,
    const dictionary& dict,
    const fvMesh& mesh
)
:
    AveragingMethod<Type>(io, dict, mesh, size(mesh)),
    volumeCell_(mesh.V()),
    cellValue_
    (
        IOobject
        (
            this->name() + ":cellValue",
            this->time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensioned<Type>("zero", dimless, pTraits<Type>::zero)
    ),
    cellCounter_
    (
        IOobject
        (
            this->name() + ":cellCounter",
            this->time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar("zero", dimless, 0.0)
    ),
    pointVolume_
    (
        IOobject
        (
            this->name() + ":pointVolume",
            this->time().timeName(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        pointMesh(mesh),
        dimensionedScalar("zero", dimless, 0.0)
    ),
    cellUniqueList_(mesh.nCells()),
    cellTotalSize_(mesh.nCells(),0),
    S_(vector::zero),
    pL_(mesh.nCells()),
    pOrder_(mesh.nCells())
{
    //initialise coefficients lists
    forAll(mesh.cells(),cI)
    {
        pL_[cI].setSize(8, vector::zero);
        pOrder_[cI].setSize(8,0);
    }

    //- order quadralaterial points in correct order
    order_hexa();

    //- calculate point volume weights to correctly account for boundary
    //  corner nodes
    forAll(mesh.cells(),cI)
    {
        labelList cellLinear;
        labelList pID = pOrder_[cI];

        forAll (pID,p)
        {
            const labelList& cID = this->mesh_.pointCells(pID[p]);
            cellLinear.append(cID);

            pointVolume_[pID[p]] += 0.5*0.5*0.5;
        }

        Foam::sort(cellLinear);
        labelList cellUnique;
        Foam::uniqueOrder(cellLinear,cellUnique);

        cellTotalSize_[cI] = cellUnique.size();

        cellUniqueList_[cI].setSize(cellTotalSize_[cI]);

        forAll(cellUnique,I)
        {
            cellUniqueList_[cI][I] = cellLinear[cellUnique[I]];
        }
    }
}

template<class Type>
Foam::AveragingMethods::Linear<Type>::Linear
(
    const Linear<Type>& am
)
:
    AveragingMethod<Type>(am),
    volumeCell_(am.volumeCell_),
    cellValue_(am.cellValue_),
    cellCounter_(am.cellCounter_),
    pointVolume_(am.pointVolume_),
    cellUniqueList_(am.cellUniqueList_),
    cellTotalSize_(am.cellTotalSize_),
    S_(am.S_),
    pL_(am.pL_),
    pOrder_(am.pOrder_)
{

}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class Type>
Foam::AveragingMethods::Linear<Type>::~Linear()
{}


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Linear<Type>::order_hexa
(
)
{
    //- loops through cells, finds minimum y face, calculates opposite face
    //  and builds point list in order required by XtoL

    //- loop through cells
    forAll(this->mesh_.C(),cI)
    {
        const cell& cells = this->mesh_.cells()[cI];

        label lF      = 0;
        scalar lFz  = VGREAT;

        //- loop through cell faces
        forAll(cells, fI )
        {
            const face& faces = this->mesh_.faces()[cells[fI]];
            pointField points = faces.points(this->mesh_.points());

            scalar tmpf = 0.0;

            //- determine lowest face in y direction
            forAll(faces,pI)
            {
                tmpf     += points[pI].y();
            }

            if(tmpf < lFz)
            {
                lF  = fI;
                lFz = tmpf;
            }
        }

        const face& lFface = this->mesh_.faces()[cells[lF]];
        oppositeFace opFface = cells.opposingFace(cells[lF],this->mesh_.faces());

        pointField lFPoints = lFface.points(this->mesh_.points());
        pointField opFPoints = opFface.points(this->mesh_.points());

        forAll(lFface,pI)
        {
            pL_[cI][pI] = lFPoints[pI];
            pOrder_[cI][pI] = lFface[pI];
        }

        forAll(opFface,pI)
        {
            pL_[cI][pI+4] = opFPoints[pI];
            pOrder_[cI][pI+4] = opFface[pI];
        }

        /*Info << "Cell[" << cI << "] "
             << "points:" << pL_[cI]
             << nl
             << "pointOrder: " << pOrder_[cI]
             << endl;

        Info    << "Cell[" << cI
                    << "] Face[" << cells[lF]
                    << nl
                    << "    points:" << lFface.points(this->mesh_.points())
                    << nl
                    << "    faceNormal" << lFface.normal(this->mesh_.points())/lFface.mag(this->mesh_.points())
                    << endl;

        oppositeFace opFface = cells.opposingFace(cells[lF],this->mesh_.faces());

        Info    << "Opposite Face "
                    << nl
                    << "    points:" << opFface.points(this->mesh_.points())
                    << nl
                    << "    faceNormal" << opFface.normal(this->mesh_.points())/opFface.mag(this->mesh_.points())
                    << endl;*/


    }
}

template<class Type>
void Foam::AveragingMethods::Linear<Type>::XtoL
(
    const point position,
    const label& cI
) const
{
    //computes the inverse trilinear map from [0,1]^3 to the arbritrary
    //quadralateral defined by the points p0-7, where points are ordered
    //consistent with
    //p0~(0,0,0), p1~(1,0,0), p2~(1,0,1), p3~(0,0,1),
    //p4~(0,1,0), p5~(1,1,0), p6~(1,1,1), p7~(0,1,1),
    //Uses Gauss-Newton method.

    label iter = 30;                  //- max iteration counter
    scalar tol = 1e-10;               //- tolerance
    vector ss(0.5, 0.5, 0.5);         //- initial guess
    //- logical co-ords
    scalar dl = 0;
    scalar dm = 0;
    scalar dp = 0;

    scalar normR = 1;
    //- operators
    vector r  = vector::zero;
    vector xr = vector::zero;

    vector Js = vector::zero;
    vector Jt = vector::zero;
    vector Jw = vector::zero;

    tensor J  = tensor::zero;
    tensor JT  = tensor::zero;
    tensor A    = tensor::zero;
    tensor invA = tensor::zero;

    for(label k = 1; k<iter; k++)
    {
        dl = ss.x();
        dm = ss.y();
        dp = ss.z();

        r = pL_[cI][0]*(1-dl)*(1-dm)*(1-dp) + pL_[cI][1]*  dl * (1-dm)*(1-dp) +
            pL_[cI][2]*   dl *(1-dm)*   dp  + pL_[cI][3]*(1-dl)*(1-dm)*   dp +
            pL_[cI][4]*(1-dl)*   dm *(1-dp)  + pL_[cI][5]*   dl *   dm *(1-dp) +
            pL_[cI][6]*   dl *   dm *   dp   + pL_[cI][7]*(1-dl)*   dm*    dp - position;

        normR = sqrt(r.x()*r.x() + r.y()*r.y() + r.z()*r.z());

        if (normR < tol)
        {
            break;
        }

        // Construct Jacobian matrix
        Js = -pL_[cI][0]*(1-dm)*(1-dp) + pL_[cI][1]*(1-dm)*(1-dp) +
              pL_[cI][2]*(1-dm)*   dp  - pL_[cI][3]*(1-dm)*   dp  +
             -pL_[cI][4]*   dm*(1-dp)  + pL_[cI][5]*   dm *(1-dp) +
              pL_[cI][6]*   dm*   dp   - pL_[cI][7]*   dm *   dp;

        Jt = -pL_[cI][0]*(1-dl)*(1-dp) - pL_[cI][1]*   dl *(1-dp) +
             -pL_[cI][2]*   dl *   dp  - pL_[cI][3]*(1-dl)*   dp  +
              pL_[cI][4]*(1-dl)*(1-dp) + pL_[cI][5]*   dl *(1-dp) +
              pL_[cI][6]*   dl *   dp  + pL_[cI][7]*(1-dl)*   dp;

        Jw = -pL_[cI][0]*(1-dl)*(1-dm) - pL_[cI][1]*  dl *(1-dm)  +
              pL_[cI][2]*   dl *(1-dm) + pL_[cI][3]*(1-dl)*(1-dm) +
             -pL_[cI][4]*(1-dl)*   dm  - pL_[cI][5]*   dl *   dm  +
              pL_[cI][6]*   dl *   dm  + pL_[cI][7]*(1-dl)*   dm ;

        J = tensor(Js.x(),Jt.x(),Jw.x(),Js.y(),Jt.y(),Jw.y(),Js.z(),Jt.z(),Jw.z());
        JT = J.T();

        A = JT & J;
        invA = inv(A);

        xr = JT & r;
        ss = ss - (invA & xr);
    }
    S_ = ss;
}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Type>
void Foam::AveragingMethods::Linear<Type>::resetFields
(
)
{
    forAll(cellValue_,cI)
    {
        cellValue_[cI] =  pTraits<Type>::zero;
    }
}

template<class Type>
void Foam::AveragingMethods::Linear<Type>::add
(
    const point position,
    const tetIndices& tetIs,
    const Type& value
)
{
    //Input is particle data, output is volume field

    const label& cI = tetIs.cell();
    //convert parcel physical co-ords to logical
    XtoL(position,cI);

    scalar dl = S_.x();
    scalar dm = S_.y();
    scalar dp = S_.z();

    //initialise useful data
    labelList pID = pOrder_[cI];
    scalar vol = volumeCell_[cI];

    List<Type> pointValue;
    pointValue.setSize(8,pTraits<Type>::zero);

    //Linear assignment
    pointValue[0] += value/(pointVolume_[pID[0]]*8*vol)*(1-dl)*(1-dm)*(1-dp);
    pointValue[1] += value/(pointVolume_[pID[1]]*8*vol)*   dl *(1-dm)*(1-dp);
    pointValue[2] += value/(pointVolume_[pID[2]]*8*vol)*   dl *(1-dm)*   dp;
    pointValue[3] += value/(pointVolume_[pID[3]]*8*vol)*(1-dl)*(1-dm)*   dp;
    pointValue[4] += value/(pointVolume_[pID[4]]*8*vol)*(1-dl)*   dm* (1-dp);
    pointValue[5] += value/(pointVolume_[pID[5]]*8*vol)*   dl *   dm* (1-dp);
    pointValue[6] += value/(pointVolume_[pID[6]]*8*vol)*   dl*    dm*    dp;
    pointValue[7] += value/(pointVolume_[pID[7]]*8*vol)*(1-dl)*   dm*    dp;

    forAll(pID,p)
    {
        //interpolate value to cell nodes
        const labelList& cID = this->mesh_.pointCells(pID[p]);

        forAll(cID,c)
        {
            //Corrects for non-uniform volumes
            cellValue_[cID[c]] += pointValue[p]*vol/(volumeCell_[cID[c]]);
        }
    }
}

template<class Type>
Foam::List<List<scalar> > Foam::AveragingMethods::Linear<Type>::stencil
(
    const point position,
    const tetIndices& tetIs
)
{

    const label& cI = tetIs.cell();

    //- calculate stencil for particle/cell combination
    XtoL(position,cI);

    scalar dl = S_.x();
    scalar dm = S_.y();
    scalar dp = S_.z();

    //initialise useful data
    labelList pID = pOrder_[cI];

    List<scalar> pointValue;

    scalar value = 1.0;
    pointValue.setSize(8,0.0);

    //Linear assignment
    pointValue[0] += value/(pointVolume_[pID[0]]*8)*(1-dl)*(1-dm)*(1-dp);
    pointValue[1] += value/(pointVolume_[pID[1]]*8)*   dl *(1-dm)*(1-dp);
    pointValue[2] += value/(pointVolume_[pID[2]]*8)*   dl *(1-dm)*   dp;
    pointValue[3] += value/(pointVolume_[pID[3]]*8)*(1-dl)*(1-dm)*   dp;
    pointValue[4] += value/(pointVolume_[pID[4]]*8)*(1-dl)*   dm* (1-dp);
    pointValue[5] += value/(pointVolume_[pID[5]]*8)*   dl *   dm* (1-dp);
    pointValue[6] += value/(pointVolume_[pID[6]]*8)*   dl*    dm*    dp;
    pointValue[7] += value/(pointVolume_[pID[7]]*8)*(1-dl)*   dm*    dp;

    List<scalar> cellValueStencil;
    cellValueStencil.setSize(cellTotalSize_[cI],0.0);
    scalar totVol = 0.0;
    forAll(cellUniqueList_[cI],cellI)
    {
        totVol += volumeCell_[cellUniqueList_[cI][cellI]];
    }

    forAll(pID,p)
    {
        //interpolate value to cell nodes
        const labelList& cID = this->mesh_.pointCells(pID[p]);

        //- add contribution to cell
        forAll(cID,c)
        {
            //Corrects for non-uniform volumes
            label index = Foam::findIndex(cellUniqueList_[cI],cID[c]);

            if(index != -1)
            {
                //cellValueStencil[index] += pointValue[p];
                cellValueStencil[index] += pointValue[p];
            }
        }
    }

    //setup result list
    List<List<scalar> > result(2);
    forAll(result,I)
    {
        result[I].setSize(cellTotalSize_[cI]);
    }

    //grab result
    forAll(cellUniqueList_[cI],I)
    {
        result[0][I] = cellUniqueList_[cI][I];
        result[1][I] = cellValueStencil[I];
    }

    return result;
}

template<class Type>
Type Foam::AveragingMethods::Linear<Type>::interpolate
(
    const point position,
    const tetIndices& tetIs
)
{
    const label& cI = tetIs.cell();

    //convert parcel physical co-ords to logical
    XtoL(position,cI);

    scalar dl = S_.x();
    scalar dm = S_.y();
    scalar dp = S_.z();

    //get ordered points for this cell
    labelList pID = pOrder_[cI];

    //reset counters
    List<Type> pointValue;
    pointValue.setSize(8,pTraits<Type>::zero);

    List<scalar> pointWeight;
    pointWeight.setSize(8,0);

    forAll(pID,p)
    {
        const labelList& cID = this->mesh_.pointCells(pID[p]);
        forAll(cID,c)
        {
            cellCounter_[cID[c]] = 0.0;
        }
    }

    //determine which cells had volumes added to them from this cell and how many per cell
    forAll(pID,p)
    {
        const labelList& cID = this->mesh_.pointCells(pID[p]);
        forAll(cID,c)
        {
            cellCounter_[cID[c]] += 1.0;
        }
    }

    //account for weight assigned to nodes from this cell
    pointWeight[0] = (1-dl)*(1-dm)*(1-dp);
    pointWeight[1] =    dl *(1-dm)*(1-dp);
    pointWeight[2] =    dl *(1-dm)*   dp;
    pointWeight[3] = (1-dl)*(1-dm)*   dp;
    pointWeight[4] = (1-dl)*   dm* (1-dp);
    pointWeight[5] =    dl *   dm* (1-dp);
    pointWeight[6] =    dl*    dm*    dp;
    pointWeight[7] = (1-dl)*   dm*    dp;

    scalar vol = volumeCell_[cI];

    forAll(pID,p)
    {
        //interpolate value of cells surrounding point node to node
        const labelList& cID = this->mesh_.pointCells(pID[p]);

        forAll(cID,c)
        {
            //Accounts for variable volumes
            if(cellCounter_[cID[c]] > VSMALL && pointWeight[p] > VSMALL)
            {
                pointValue[p] += cellValue_[cID[c]]/pointWeight[p]/cellCounter_[cID[c]] * volumeCell_[cID[c]]/vol;
            }
        }
    }

    //Linear charge assignment using corrected pointValues
    return
          (pointValue[0]*(1-dl)*(1-dm)*(1-dp)
        + pointValue[1]*   dl *(1-dm)*(1-dp)
        + pointValue[2]*   dl *(1-dm)*   dp
        + pointValue[3]*(1-dl)*(1-dm)*   dp
        + pointValue[4]*(1-dl)*   dm*(1-dp)
        + pointValue[5]*   dl *   dm*(1-dp)
        + pointValue[6]*   dl *   dm*   dp
        + pointValue[7]*(1-dl)*   dm*   dp)/cellTotalSize_[cI];

    /*
    forAll(pID,p)
    {
        //interpolate value of cells surrounding point node to node
        const labelList& cID = this->mesh_.pointCells(pID[p]);

        forAll(cID,c)
        {
            //Accounts for variable volumes
            pointValue[p] += cellValue_[cID[c]]/(cellCounter_[cID[c]]);
        }

        //Account for weights
        pointValue[p] /= pointWeight[p];

    }

    //Linear charge assignment using corrected pointValues
    return
          (pointValue[0]*(1-dl)*(1-dm)*(1-dp)
        + pointValue[1]*   dl *(1-dm)*(1-dp)
        + pointValue[2]*   dl *(1-dm)*   dp
        + pointValue[3]*(1-dl)*(1-dm)*   dp
        + pointValue[4]*(1-dl)*   dm*(1-dp)
        + pointValue[5]*   dl *   dm*(1-dp)
        + pointValue[6]*   dl *   dm*   dp
        + pointValue[7]*(1-dl)*   dm*   dp)/cellTotalSize_[cI];
    */
}


template<class Type>
typename Foam::AveragingMethods::Linear<Type>::TypeGrad
Foam::AveragingMethods::Linear<Type>::interpolateGrad
(
    const point position,
    const tetIndices& tetIs
) const
{
    TypeGrad temp(
                    pTraits<Type>::zero,
                    pTraits<Type>::zero,
                    pTraits<Type>::zero
                 );

    return temp;
}


template<class Type>
void Foam::AveragingMethods::Linear<Type>::average()
{
}


template<class Type>
void Foam::AveragingMethods::Linear<Type>::average
(
    const AveragingMethod<scalar>& weight
)
{
}

template<class Type>
Foam::tmp<Foam::Field<Type> >
Foam::AveragingMethods::Linear<Type>::internalField() const
{
    return tmp<Field<Type> >(cellValue_.internalField());
}

template<class Type>
void Foam::AveragingMethods::Linear<Type>::updateField
(
    GeometricField<Type, fvPatchField, volMesh> cField
)
{
    forAll(cellValue_,cI)
    {
        cellValue_[cI ]= cField[cI];
    }

}


// ************************************************************************* //
