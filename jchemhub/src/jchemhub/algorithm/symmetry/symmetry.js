/*
 * Copyright 2010 Tim Vandermeersch
 *
 * Licensed under the Apache License, Version 2.0 (the "License"); 
 * you may not use this file except in compliance with the License. 
 * You may obtain a copy of the License at http://www.apache.org/licenses/LICENSE-2.0 
 * 
 * Unless required by applicable law or agreed to in writing, software 
 * distributed under the License is distributed on an "AS IS" BASIS, 
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. 
 * See the License for the specific language governing permissions and 
 * limitations under the License.
 */

goog.provide('jchemhub.symmetry.Symmetry');

goog.require('goog.array');

(function() {
 
    /**
     * Create an initial distance matrix. This is  the D matrix from the paper.
     */
    function createWeightMatrix(molecule, n) {
        var matrix = [];
        for (var i = 0; i < n; i++) {
            var row = [];
            for (var j = 0; j < n; j++) {
                if (i === j) {
                    row.push(0);
                } else if (molecule.findBond(molecule.getAtom(i), molecule.getAtom(j))) {
                    row.push(1);
                } else {
                    row.push(99999999);
                }
            }
            matrix.push(row);
        }
        return matrix;
    }


    /**
     * Floyd-Warshall algorithm for distance matrix.
     */
    function distanceMatrix(molecule) {
        var n = molecule.countAtoms();
        var D = createWeightMatrix(molecule, n);
        for (var k = 0; k < n; k++) {
            for (var i = 0; i < n; i++) {
                for (var j = 0; j < n; j++) {
                    D[i][j] = Math.min( D[i][j], D[i][k] + D[k][j] );
                }
            }
        }
        return D;    
    }

    /**
     * Graph Theoretical Distance. This is the distance from an atom to all other
     * atoms in the molecule. 
     */
    function graphTheoreticalDistance(molecule) {
        var n = molecule.countAtoms();
        var D = distanceMatrix(molecule);
        var GTD = [];
        for (var i = 0; i < n; i++) {
            var distance = 0;
            for (var j = 0; j < n; j++) {
                distance += D[i][j];            
            }
            GTD.push(distance);
        }
        return GTD;
    }

    function graphInvariants(molecule) {
        var n = molecule.countAtoms();
        var GTD = graphTheoreticalDistance(molecule);
        var invariants = [];
        for (var i = 0; i < n; i++) {
            var invariant = GTD[i]; // Add more invariants: atomic number, heavy valence, ring, charge 
            invariants.push(invariant);
        }
        return invariants;
    }

    function extendSymmetryClasses(molecule, symmetry_classes) {
        var n = molecule.countAtoms();
        var newClasses = [];
        for (var i = 0; i < n; i++) {
            var nbrs = molecule.getAtom(i).getNeighbors();
            var sum = 0;
            for (var j = 0, jl = nbrs.length; j < jl; j++) {
                var nbrIndex = molecule.indexOfAtom(nbrs[j]);
                sum += symmetry_classes[nbrIndex];
            }
            newClasses.push(sum);
        }
        symmetry_classes = newClasses;    
    }

    jchemhub.symmetry.countSymmetryClasses = function(symmetry_classes) {
        var n = symmetry_classes.length;
        var unique = [];
        for (var i = 0; i < n; i++) {
            if (goog.array.indexOf(unique, symmetry_classes[i]) === -1) {
                unique.push(symmetry_classes[i]);
            } 
        }
        return unique.length;
    };

    function renumberSymmetryClasses(symmetry_classes) {
        var n = symmetry_classes.length;
        var unique = [];
        for (var i = 0; i < n; i++) {
            if (goog.array.indexOf(unique, symmetry_classes[i]) === -1) {
                unique.push(symmetry_classes[i]);
            } 
        }
        unique.sort(function(a,b){return a-b;});
        var newClasses = [];
        for (i = 0; i < n; i++) {
            newClasses.push(goog.array.indexOf(unique, symmetry_classes[i])+1);
        }
        return newClasses;    
    }

    /**
     * @namespace jchemhub
     * @name jchemhub
     */
    /**
     * @namespace jchemhub.symmetry
     * @name jchemhub.symmetry
     */

    /**
     * Class for finding symmetry classes for molecules.
     * @class
     */
    jchemhub.symmetry.Symmetry = {};

    /**
     * Find the symmetry classes for a molecule.
     * @param molecule The input molecule.
     * @return {Array.<number>} The symmetry classes sorted by atom index.
     */
    jchemhub.symmetry.Symmetry.findSymmetryClasses = function(/**jchemhub.model.Molecule*/molecule) {
        var symmetry_classes = graphInvariants(molecule);

        var lastCount = jchemhub.symmetry.countSymmetryClasses(symmetry_classes);
        while (true) {

            extendSymmetryClasses(molecule, symmetry_classes);
            
            var newCount = jchemhub.symmetry.countSymmetryClasses(symmetry_classes);
            if (newCount === lastCount) {
                break;
            }
        }

        symmetry_classes = renumberSymmetryClasses(symmetry_classes);

        return symmetry_classes;    
    };

 
}());
