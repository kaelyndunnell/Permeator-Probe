/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

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

#include "codedFixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "read.H"

//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = b1a6224189253c5713b5d843663ab0769b73c881
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void parabolicVelocity_b1a6224189253c5713b5d843663ab0769b73c881(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    parabolicVelocityFixedValueFvPatchVectorField
);


const char* const parabolicVelocityFixedValueFvPatchVectorField::SHA1sum =
    "b1a6224189253c5713b5d843663ab0769b73c881";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct parabolicVelocity sha1: b1a6224189253c5713b5d843663ab0769b73c881"
            " from patch/dictionary\n";
    }
}


parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const parabolicVelocityFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct parabolicVelocity sha1: b1a6224189253c5713b5d843663ab0769b73c881"
            " from patch/DimensionedField/mapper\n";
    }
}


parabolicVelocityFixedValueFvPatchVectorField::
parabolicVelocityFixedValueFvPatchVectorField
(
    const parabolicVelocityFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct parabolicVelocity sha1: b1a6224189253c5713b5d843663ab0769b73c881 "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

parabolicVelocityFixedValueFvPatchVectorField::
~parabolicVelocityFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy parabolicVelocity sha1: b1a6224189253c5713b5d843663ab0769b73c881\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void parabolicVelocityFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs parabolicVelocity sha1: b1a6224189253c5713b5d843663ab0769b73c881\n";
    }

//{{{ begin code
    #line 44 "/home/ubuntu/OpenFOAM/kae-13/permeator-probe/0/U/inlet"

            const vectorField& Cf = patch().Cf();
            vectorField& field = *this;

            const scalar c = 0; 
            const scalar r = 0.065;
            const scalar Umax = 0.0076;


            forAll(Cf, faceI)
            {
                const scalar x = Cf[faceI][0];
                const scalar z = Cf[faceI][2];  // flow is in y, so x-z is cross-section
                const scalar rCoord = sqrt(pow(x - c, 2) + pow(z - c, 2));

                field[faceI] = vector(0, Umax * (1 - 0.5 * pow(rCoord / r, 2)), 0);



            }
        
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

