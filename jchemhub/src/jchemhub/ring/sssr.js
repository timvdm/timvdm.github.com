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

goog.provide('jchemhub.ring.SSSRFinder');

goog.require('goog.structs.Set');
goog.require('goog.array');
goog.require('jchemhub.ring.Ring');

/*
function debug(text, noNewLine) {
    var logDiv = document.getElementById("logDiv");
    if (!logDiv) {
        var body = document.getElementsByTagName("body")[0];
        logDiv = document.createElement("div");
        logDiv.id = "logDiv";
        body.appendChild(logDiv);
    }
    if (noNewLine) {
        logDiv.innerHTML += text;
    } else {
        logDiv.innerHTML += text + "<br>";
    }
};
*/

(function() {

    /**
     * Smallest Set of Smallest rings.
     *
     * A robust method for searching the smallest set of smallest rings with a 
     * path-included distance matrix, Chang Lee, PNAS, Oct. 13, 2009, vol. 16,
     * no. 41, pages 17355-17358
     *
     * http://www.pnas.org/content/106/41/17355.full
     * http://www.pnas.org/content/106/41/17355/suppl/DCSupplemental
     */


    /**
     * Make a deep copy of an array
     */
    function deepCopy(arr) {
        var newArray = [];
        for (var i = 0, li = arr.length; i < li; i++) {
            var item = arr[i];
            if (item instanceof Array) {
                newArray.push(deepCopy(item));
            } else {
                newArray.push(item);
            }
        }
        return newArray;
    }

    /**
     * Debug helper function
     */
    function matrixToHTML(matrix) {
        var text = "";
        var n = matrix.length;
        for (var i = 0; i < n; i++) {
            for (var j = 0; j < n; j++) {
                text += JSON.stringify(matrix[i][j]) + " ";
            }
            text += "<br>";
        }
        return text; 
    }

    /**
     * Create a n x n matrix with all elements set to 0.
     */
    function createEmptyMatrix(n) {
        var matrix = [];
        for (var i = 0; i < n; i++) {
            var row = [];
            for (var j = 0; j < n; j++) {
                row.push(0);
            }
            matrix.push(row);
        }
        return matrix;
    }

    /**
     * Create an initial distance matrix. This is  the D matrix from the paper.
     */
    function createWeightMatrix(molecule, n) {
        var matrix = [];
        for (var i = 0; i < n; i++) {
            var row = [];
            for (var j = 0; j < n; j++) {
                if (i == j) {
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
     * Create an empty Path-Included Distance matrix. This is an n x n matrix 
     * with all elements set to an empty list (i.e. there is no path between 
     * i and j)
     */
    function createEmptyPIDMatrix(n) {
        var matrix = [];
        for (var i = 0; i < n; i++) {
            var row = [];
            for (var j = 0; j < n; j++) {
                row.push([]);
            }
            matrix.push(row);
        }
        return matrix;

    }

    /**
     * Create the initial Pe matrix. The supplementary information 
     * indicates this should be an empty PID matrix, this makes no sense though.
     * We initialize it with adding all bonds as paths of length 1 when atom i 
     * and j are connected.
     */
    function createPIDMatrix(molecule, n) {
        var matrix = [];
        for (var i = 0; i < n; i++) {
            var row = [];
            for (var j = 0; j < n; j++) {
                var beginAtom = molecule.getAtom(i);
                var endAtom = molecule.getAtom(j);
                var bond = molecule.findBond(beginAtom, endAtom);
                if (bond) {
                    row.push([[molecule.indexOfBond(bond)]]);
                } else {
                    row.push([]);
                }
            }
            matrix.push(row);
        }
        return matrix;

    }

    /**
     * Update a path included distance matrix element (i.e. lhs) by merging 
     * p1 with p2. 
     */
    function appendPath(lhs, p1, p2) {
        if (!lhs.length) {
            lhs[0] = p1[0].concat(p2[0]);
        } else {
            lhs.push(p1[0].concat(p2[0]));
        }
    }


    /**
     * Create the two Path-Included Distance matrices (Pe and Pe') and the 
     * distance matrix D. This is Algorithm 1 in the supplementary information
     * that comes with the paper.
     */
    function makePIDMatrixes(molecule) {
        var n = molecule.countAtoms();
        var D = createWeightMatrix(molecule, n);
        var Pe1 = createPIDMatrix(molecule, n); // Pe
        var Pe2 = createEmptyPIDMatrix(n); // Pe'
        var lastD = D;

        //debug("Pe =<br>" + matrixToHTML(Pe1));
        //debug("Pe' =<br>" + matrixToHTML(Pe2));

        for (var k = 0; k < n; k++) {
            for (var i = 0; i < n; i++) {
                for (var j = 0; j < n; j++) {
                    var lastPathLength = lastD[i][j];
                    var pathLength = lastD[i][k] + lastD[k][j];
                    var path1 = Pe1[i][k];
                    var path2 = Pe1[k][j];

                    // ignore invalid paths
                    if (pathLength == 100000000) { continue; }

                    if (lastPathLength > pathLength) {
                        if (lastPathLength == pathLength + 1) {
                            // a new shortest path = previous shortest path -1 => Pe' <- Pe
                            Pe2[i][j] = deepCopy(Pe1[i][j]);
                        } else {
                            Pe2[i][j] = [];
                        }

                        // a new shortest path is found
                        D[i][j] = pathLength;
                        Pe1[i][j] = [ path1[0].concat(path2[0]) ]; // change path
                    } else if (lastPathLength == pathLength) {
                        // another shortest path is found
                        if (path1.length && path2.length) {
                            appendPath(Pe1[i][j], path1, path2); // append path
                        }
                    } else if (lastPathLength == pathLength - 1) {
                        // shortest path + 1 found
                        appendPath(Pe2[i][j], path1, path2); // append path
                    } else {
                        D[i][j] = lastD[i][j];
                    }
                }
            }
        }

        //debug("D =<br>" + matrixToHTML(D));
        //debug("Pe =<br>" + matrixToHTML(Pe1));
        //debug("Pe' =<br>" + matrixToHTML(Pe2));

        return { "D": D, "Pe1": Pe1, "Pe2": Pe2 };
    };

    /**
     * Sort function to sort the set of candidates by increasing Cnum (i.e. ring size).
     */
    function sortByCnum(a, b) {
        return a.Cnum - b.Cnum;
    }
  
    /**
     * Compute the set of ring candidates using the distance matrix and the
     * two path-included distance matrices. This is algorithm 2 in supplementary
     * information.
     */
    function makeCandidateSet(D, Pe1, Pe2) {
        var n = D.length;
        var Cset = [];

        for (var i = 0; i < n; i++) {
            for (var j = 0; j < n; j++) {
                if (D[i][j] == 0 || (Pe1[i][j].length == 1 && Pe2[i][j].length == 0)) {
                    continue; // skip degenerate cases
                } else {
                    var Cnum;
                    if (Pe2[i][j].length) {
                        Cnum = 2 * (D[i][j] + 0.5); // odd ring candidate        
                    } else {
                        Cnum = 2 * D[i][j]; // even ring candidate
                    }
                    Cset.push({ "Cnum": Cnum, "Pe1": Pe1[i][j], "Pe2": Pe2[i][j] });
                }
            }    
        }

        // sort the candidates by increasing ring size
        Cset.sort(sortByCnum);

        //for (var i = 0; i < Cset.length; i++) {
        //  var C = Cset[i];
        //  debug("Cset <= " + C.Cnum + " + " + JSON.stringify(C.Pe1) + " + " + JSON.stringify(C.Pe2));
        //}

        return Cset;
    };


    /**
     * Check if a candidate is already part of the SSSR set. A ring is considered
     * to be in the set if the set contains an identical ring or a ring that is a 
     * subset of the candidate. This is the XOR function from the paper.
     */ 
    function isCandidateInSet(C, Csssr) {
        for (var i = 0, li = Csssr.length; i < li; i++) {
            var sssr = Csssr[i];
            if (C.length < sssr) {
                continue;
            }
            var candidateContainsRing = true;
            for (var j = 0, lj = sssr.length; j < lj; j++) {
                if (goog.array.indexOf(C, sssr[j]) == -1) {
                    candidateContainsRing = false;
                }
            }
            if (candidateContainsRing)
                return true;
        }
        return false;
    }


    /**
     * Search the candidates to find the Smallest Set of Smallest rings. This
     * is algorithm 3 from the supplementary information.
     */
    function candidateSearch(Cset, nsssr) {
        var ringIndex = 0;
        var Csssr = [];
        for (var i = 0, li = Cset.length; i < li; i++) {
            //debug("Cset <= " + Cset[i].Cnum + " + " + JSON.stringify(Cset[i].Pe1) + " + " + JSON.stringify(Cset[i].Pe2));
            if (Cset[i].Cnum % 2) {
                // odd ring
                for (var j = 0, lj = Cset[i].Pe2.length; j < lj; j++) {
                    var C = Cset[i].Pe1[0].concat(Cset[i].Pe2[j]);
                    //debug("C = " + JSON.stringify(C));
                    if (!isCandidateInSet(C, Csssr)) {
                        Csssr.push(C);
                        ringIndex++;
                    }
                    if (ringIndex == nsssr) {
                        return Csssr;
                    }
                }
            } else {
                // even ring
                for (var j = 0, lj = Cset[i].Pe1.length - 1; j < lj; j++) {
                    var C = Cset[i].Pe1[j].concat(Cset[i].Pe1[j+1]);
                    //debug("C = " + JSON.stringify(C));
                    if (!isCandidateInSet(C, Csssr)) {
                        Csssr.push(C);
                        ringIndex++;
                    }
                    if (ringIndex == nsssr) {
                        return Csssr;
                    }
                }
            }
        }

        return Csssr;
    };

    /**
     * Find the Smallest Set of Smallest rings.
     */
    jchemhub.ring.findSSSR = function(molecule) {
        var matrices = makePIDMatrixes(molecule)
        var Cset = makeCandidateSet(matrices.D, matrices.Pe1, matrices.Pe2);
        var nsssr = molecule.countBonds() - molecule.countAtoms() + 1;
        var indexes = candidateSearch(Cset, nsssr);

        var rings = [];
        for (var i = 0, li = indexes.length; i < li; i++) {
            var index = indexes[i];
            var atoms = [];
            var bonds = [];
            for (var j = 0, lj = index.length; j < lj; j++) {
                var bond = molecule.getBond(index[j]);
                bonds.push(bond);
                if (goog.array.indexOf(atoms, bond.source) == -1) {
                    atoms.push(bond.source);
                }
                if (goog.array.indexOf(atoms, bond.target) == -1) {
                    atoms.push(bond.target);
                }
            }

            rings.push(new jchemhub.ring.Ring(atoms, bonds));
        }

        return rings;
    }

    /*
    jchemhub.ring._makePIDMatrixes = makePIDMatrixes;
    jchemhub.ring._makeCandidateSet = makeCandidateSet;
    jchemhub.ring._candidateSearch = candidateSearch;
    */

})();
