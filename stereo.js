/**
 * Copyright 2010 Tim Vandermeersch
 * 
 * Classify stereo bonds according to IUPAC recommendations 2006.
 */

"use strict";

(function() {

function vector2(x, y) 
{
  this.x = x;
  this.y = y;

  this.length = function length() {
    return Math.sqrt(this.x * this.x + this.y * this.y)
  }
  this.add = function add(other) {
    return new vector2(this.x + other.x, this.y + other.y);
  }
  this.subtract = function subtract(other) {
    return new vector2(this.x - other.x, this.y - other.y);
  }
  this.multiply = function multiply(value) {
    return new vector2(this.x * value, this.y * value);
  }
  this.dot = function dot(other) {
    return this.x * other.x + this.y * other.y;
  }
  this.cross = function cross(other) {
    return new vector2(this.x * other.y, this.y * other.x);
  }
  this.normalize = function normalize() {
    l = this.length();
    this.x /= l;
    this.y /= l;
  }
  this.normalized = function() {
    l = this.length();
    return new vector2(this.x / l, this.y / l);
  }
}
function getAngle(v1, v2) {
  var n1 = v1.normalized();
  var n2 = v2.normalized();
  var val = n1.dot(n2);
  return Math.acos(val);
}
function getAngleSign(v1, v2) {
  return v1.x * v2.y - v1.y * v2.x;
}
function getTriangleSign(a, b, c)
{
  return (a.x - c.x) * (b.y - c.y) - (a.y - c.y) * (b.x - c.x);
}
function radToDeg(radians) {
  return 180 * radians / Math.PI;
}









  window.JSChem.stereo = {};

  window.JSChem.stereo.Stereo = function(molecule, stereogenic_units) {
    
    this.molecule = molecule;
    this.units = stereogenic_units;

    this.isTetrahedral = function(atom) {
      var atomIndex = window.JSChem.getAtomIndex(this.molecule, atom);
      var numUnits = this.units.length;
      for (var i = 0; i < numUnits; i++) {
        var unit = this.units[i];
        if (unit.type == "tetrahedral" && unit.atomIndex == atomIndex) {
          return true;
        }      
      }
      return false;
    };

    this.evaluate2D = function(sssr) {

      // 
      // Helper functions
      //
      function isInRing(molecule, bond, sssr) {
        var numRings = sssr.length;
        for (var i = 0; i < numRings; i++) {
          var ring = sssr[i];
          var ringSize = ring.length;
          for (var j = 1; j < ringSize; j++) {
            var begin = molecule.atoms[ring[j-1]];
            var end = molecule.atoms[ring[j]];
            if (window.JSChem.getBondBetween(molecule, begin, end) == bond) {
              return true;            
            }
          }
          var begin = molecule.atoms[ring[0]];
          var end = molecule.atoms[ring[ringSize-1]];
          if (window.JSChem.getBondBetween(molecule, begin, end) == bond) {
            return true;            
          }
        }
        return false;
      }
      function computeIndexes(sortedAngles, index1, index2) {
        var i1, i2;
        var numAngles = sortedAngles.length;
        for (var i = 0; i < numAngles; i++) {
          var angle = sortedAngles[i];
          if (angle.index == index1) {
            i1 = i;
          } else if (angle.index == index2) {
            i2 = i;
          }
        }
        return [i1, i2];
      }
      function computeAngle(sortedAngles, index1, index2) {
        var i = computeIndexes(sortedAngles, index1, index2);
        var imax = Math.max(i[0], i[1]);
        var imin = Math.min(i[0], i[1]);
        var maxAngle = sortedAngles[imax].angle;
        var minAngle = sortedAngles[imin].angle;
        var angle = maxAngle - minAngle;
        if (angle > 180) {
          angle = 360 - angle;
        }
        return angle;      
      }
      function isAdjacent(sortedAngles, index1, index2) {
        var i = computeIndexes(sortedAngles, index1, index2);
        var max = Math.max(i[0], i[1]);
        var min = Math.min(i[0], i[1]);
        if (max - min == 1) {
          return true;
        }
        if (min == 0 && max == sortedAngles.length - 1) {
          return true;
        }
        return false;
      }
      function matchAngle(ref, angle) {
        if (angle > ref - 10.0 && angle < ref + 10.0) {
          return true;
        }
        return false;      
      }
      function computeGap(sortedAngles) {
        var gaps = [];
        var numAngles = sortedAngles.length;
        for (var i = 1; i < numAngles; i++) {
          gaps.push( sortedAngles[i].angle - sortedAngles[i-1].angle );          
        }
        gaps.push( 360.0 - sortedAngles[numAngles-1].angle + sortedAngles[0].angle );
        function sortDescending(a, b) {
          return b - a;
        }
        gaps.sort(sortDescending);
        return gaps[0];
      }
      function isSpiroAtom(atomIndex, sssr) {
        // todo: performance... check if atom is in any ring first...
        var numRings = sssr.length;
        for (var i = 0; i < numRings; i++) {
          var iRing = sssr[i];
          var iSize = iRing.length;
          for (var j = 0; j < numRings; j++) {
            if (i == j) {
              continue;
            }              
            var jRing = sssr[j];
            var jSize = jRing.length;
            
            var intersection = [];
            for (var k = 0; k < iSize; k++) {
              for (var l = 0; l < jSize; l++) {
                if (iRing[k] == jRing[l]) {
                  intersection.push(iRing[k]);
                }
              }            
            }

            if (intersection.length == 1) {
              if (intersection[0] == atomIndex) {
                return true;
              }
            }
          }        
        }
      }
      

      var result = [];
      //
      // Check for stereo bonds between two tetrahedral centers
      // - this is only valid when the 4 neighbor atom are also tetrahedral stereo centers
      //
      // Also check for inverted (a.k.a. perspective) hash bonds
      //
      var perspectiveHash = false;
      var invalidBonds = [];
      var numBonds = this.molecule.bonds.length;
      for (var i = 0; i < numBonds; i++) {
        var bond = this.molecule.bonds[i];
        if (!bond.type.length) {
          continue; // consider only wedge/hash bonds
        }
        var begin = this.molecule.atoms[bond.begin];
        var end = this.molecule.atoms[bond.end];
        // wide end of the wedge is located at the end atom
        var beginIsTetrahedral = this.isTetrahedral(begin);
        var endIsTetrahedral = this.isTetrahedral(end);
        if (beginIsTetrahedral && endIsTetrahedral) {
          // count stereogenic neighbors for begin atom
          var beginNbrs = JSChem.getAtomNeighbors(molecule, begin);
          var beginNbrsTetrahedral = 0;
          for (var j = 0, numBeginNbrs = beginNbrs.length; j < numBeginNbrs; j++) {
            if (this.isTetrahedral(beginNbrs[j])) {
              beginNbrsTetrahedral++;                
            }
          }
          // count stereogenic neighbors for end atom
          var endNbrs = JSChem.getAtomNeighbors(molecule, end);
          var endNbrsTetrahedral = 0;
          for (var j = 0, numEndNbrs = endNbrs.length; j < numEndNbrs; j++) {
            if (this.isTetrahedral(endNbrs[j])) {
              endNbrsTetrahedral++;                
            }
          }

          if (beginNbrsTetrahedral == 4 || endNbrsTetrahedral == 4) {
            // valid exception...
          } else {
            // stereo bonds between stereo centers should be avoided
            var bondIndex = JSChem.getBondIndex(this.molecule, bond);
            //JSChem.print("bond " + bondIndex + ": Not Acceptable (stereo bonds between stereo centers should be avoided, see ST-0.5)");
            result.push({ "type": "bond", "index": bondIndex, "reference": "ST-0.5", 
                "text": "Not Acceptable (stereo bonds between stereo centers should be avoided)", "ambiguous": true});
            invalidBonds.push(bond);
          }
        } else if (endIsTetrahedral) {
          if (bond.type == "hash") {
            // not really a problem yet, mixing both styles is not acceptable though
            perspectiveHash = true;
          }
        } else if (!beginIsTetrahedral) {
          // wedge/hash bond with no stereo center, probably for indicating perspective...
          var bondIndex = JSChem.getBondIndex(this.molecule, bond);
          JSChem.print("bond " + bondIndex + ": stereo bond without stereo center, to indicate perspective?");
        }
      }

      //
      // Check for stereo bonds in rings
      // - should be avoided
      // - valid:
      //   - plain bonds are to another stereo center (i.e. avoiding stereo bond between stereo centers has a higher priority)
      //   - spiro ring fusions
      // 
      var numAtoms = this.molecule.atoms.length;
      for (var i = 0; i < numAtoms; i++) {
        var atom = this.molecule.atoms[i];
        if (this.isTetrahedral(atom)) {
          if (!isSpiroAtom(JSChem.getAtomIndex(this.molecule, atom), sssr)) {
            var bonds = window.JSChem.getAtomBonds(this.molecule, atom);
            var nbrs = window.JSChem.getAtomNeighbors(this.molecule, atom);
            var nbrIsTetrahedral = [];
            var plain = [], wedge = [], hash = [];
            for (var j = 0; j < numNbrs; j++) {
              nbrIsTetrahedral.push(this.isTetrahedral(nbrs[j]));
              var bond = bonds[j];
              if (bond.type == "") {
                plain.push(j);
              } else if (bond.type == "wedge") {
                wedge.push(j);
              } else if (bond.type == "hash") {
                hash.push(j);
                if (bond.begin == JSChem.getAtomIndex(this.molecule, nbrs[j])) {
                  JSChem.print("Multiple types of hashed bonds should be avoided. ", true);
                  JSChem.print("The narrow end should either be at the tetrahedral stereo center or reverse (a.k.a. perspective hashes).");
                }
              }
            }

            if (plain.length == 2 && wedge.length == 1 && hash.length == 1) {
              var wedgeIsInRing = isInRing(this.molecule, bonds[wedge[0]], sssr);
              var hashIsInRing = isInRing(this.molecule, bonds[hash[0]], sssr);
              if (wedgeIsInRing && hashIsInRing) {
                if (!nbrIsTetrahedral[plain[0]] || !nbrIsTetrahedral[plain[1]]) {
                  // stereo bonds in ring should be avoided
                  var bondIndex = JSChem.getBondIndex(this.molecule, bond);
                  //JSChem.print("bond " + bondIndex + ": Not Acceptable (stereo bonds in rings should be avoided, see ST-0.5)");
                  result.push({ "type": "bond", "index": bondIndex, "reference": "ST-0.5", 
                      "text": "Not Acceptable (stereo bonds in rings should be avoided, see ST-0.5)", "ambiguous": true});
                  if (wedgeIsInRing) {
                    invalidBonds.push(bonds[wedge[0]]);
                  }
                  if (hashIsInRing) {
                    invalidBonds.push(bonds[hash[0]]);
                  }
                }
              }
            }
          } else {
            // ST-1.3.4
            // spiro atom
            // todo: check for distorted ring
          }
        }
      }
 
      //
      // Start the per atom tetrahedral checking...
      //
      var xAxis = new vector2(1, 0);
      var numAtoms = this.molecule.atoms.length;
      for (var i = 0; i < numAtoms; i++) {
        var atom = this.molecule.atoms[i];
        if (this.isTetrahedral(atom)) {
          var bonds = window.JSChem.getAtomBonds(this.molecule, atom);
          var nbrs = window.JSChem.getAtomNeighbors(this.molecule, atom);
          var atomPos = new vector2(atom.x, atom.y);
          var numNbrs = nbrs.length;
          // save some variables for easy access
          var skip = false;
          var bondTypes = [];
          var nbrIsInRing = [];
          var angles = [];
          var sortedAngles = [];
          var bondVectors = [];
          var numPlain = 0;
          var numWedge = 0;
          var numHash = 0;
          for (var j = 0; j < numNbrs; j++) {
            var bond = bonds[j];
            if (invalidBonds.indexOf(bond) != -1) {
              skip = true;
            }
            var nbr = nbrs[j];
            var type = bond.type.length ? bond.type : "plain" ;
            bondTypes.push( type );
            nbrIsInRing.push(isInRing(this.molecule, bond, sssr));
            var nbrPos = new vector2(nbr.x, nbr.y);
            var bv = nbrPos.subtract(atomPos);
            var angle = radToDeg(getAngle(xAxis, bv));
            if (getTriangleSign(new vector2(0.0, 0.0), xAxis, bv) > 10.0) {
              angle = 360.0 - angle;
            }
            angles.push(angle);
            sortedAngles.push({ "angle": angle, "index": j });
            bondVectors.push(bv);

            if (type == "plain") {
              numPlain++;
            } else if (type == "wedge") {
              numWedge++;
            } else if (type == "hash") {
              numHash++;
            }
          }

          if (skip) {
            continue;
          }

          function sortByAngle(a, b) {
            return a.angle - b.angle;
          }
          sortedAngles.sort(sortByAngle);

          //
          // ST-1.1 Tetrahedral configurations depicted with four explicit bonds
          //
          if (numNbrs == 4) {
            if (numPlain == 2 && numWedge == 1 && numHash == 1) {
              var plain = [];
              var wedge, hash;
              for (var j = 0; j < 4; j++) {
                if (bondTypes[j] == "plain") {
                  plain.push(j);
                } else if (bondTypes[j] == "wedge") {
                  wedge = j;
                } else if (bondTypes[j] == "hash") {
                  hash = j;              
                }
              }

              var plainAngle = computeAngle(sortedAngles, plain[0], plain[1]);
              var otherAngle = computeAngle(sortedAngles, wedge, hash);
              
              if (isAdjacent(sortedAngles, plain[0], plain[1])) {
                var largestGap = computeGap(sortedAngles);
                if (largestGap > 190.0) {
                  // ST-1.1.14
                  //wrong
                  //JSChem.print("atom " + i + ": Wrong (See ST-1.1.14)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.14", "text": "Wrong", "ambiguous": true});
                  // ambiguous!
                } else if (matchAngle(180.0, plainAngle)) {
                  // ST-1.1.13
                  // wrong
                  //JSChem.print("atom " + i + ": Wrong (See ST-1.1.13)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.13", "text": "Wrong", "ambiguous": true});
                  // ambiguous!
                } else {
                  // ST-1.1.1
                  var plainCorrect = matchAngle(120.0, plainAngle);
                  var otherCorrecr = matchAngle(60.0, otherAngle);
                  // todo: check bisectors angles
                  // todo: 1.1.1 vs 1.1.2 check fragment size
                  // preferred
                  //JSChem.print("atom " + i + ": Preferred (See ST-1.1.1)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.1", "text": "Preferred", "ambiguous": false});
                  // unambiguous -> interpret...
                }
              } else {
                if (matchAngle(180.0, plainAngle)) {
                  // ST-1.1.15
                  // wrong
                  //JSChem.print("atom " + i + ": Wrong (See ST-1.1.15)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.15", "text": "Wrong", "ambiguous": true});
                  // ambiguous!
                } else {
                  // ST-1.1.10
                  // not acceptable
                  //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.1.10)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.10", "text": "Not Acceptable", "ambiguous": true});
                  // ambiguous!
                }
              }
            } else if (numPlain == 3 && numWedge + numHash == 1) {
              var plain = [];
              var wedge_hash;
              for (var j = 0; j < 4; j++) {
                if (bondTypes[j] == "plain") {
                  plain.push(j);
                } else if (bondTypes[j] == "wedge") {
                  wedge_hash = j;
                } else if (bondTypes[j] == "hash") {
                  wedge_hash = j;              
                }
              }

              var plainAngle1 = computeAngle(sortedAngles, plain[0], plain[1]);
              var plainAngle2 = computeAngle(sortedAngles, plain[1], plain[2]);
              var plainAngle3 = computeAngle(sortedAngles, plain[0], plain[2]);

              var bothOnSameSide = 0;
              for (var j = 0; j < 3; j++) {
                var signs = [];
                for (var k = 0; k < 3; k++) {
                  if (j != k) {
                    signs.push(getTriangleSign(bondVectors[plain[j]], new vector2(0.0, 0.0), bondVectors[plain[k]]));
                  }
                }
                if (signs[0] < 10.0 && signs[1] < 10.0 || signs[0] > -10.0 && signs[1] > -10.0) {
                  bothOnSameSide++;
                }
              }

              if (plainAngle1 < 180.0 && plainAngle2 < 180.0 && plainAngle3 < 180.0 && bothOnSameSide < 2) {
                // ST-1.1.2
                matchAngle(120.0, plainAngle1);
                matchAngle(120.0, plainAngle2);
                matchAngle(120.0, plainAngle3);
                // todo: check bisect, center in bridged/fused ring system
                // todo: 1.1.1 vs 1.1.2 check fragment size
                // preferred
                //JSChem.print("atom " + i + ": Preferred (See ST-1.1.2)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.2", "text": "Preferred", "ambiguous": false});
                // unambiguous -> interpret...
              } else {
                var largestGap = computeGap(sortedAngles);
                if (largestGap > 175.0) {
                  // ST-1.1.16
                  // wrong
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.16", "text": "Wrong", "ambiguous": true});
                  // ambiguous!
                } else {
                  // todo: 2 ring bonds also acceptable...
                  if (nbrIsInRing[plain[0]] && nbrIsInRing[plain[1]] && nbrIsInRing[plain[2]]) {
                    // ST-1.1.4 exception
                    // preferred
                    //JSChem.print("atom " + i + ": Preferred (See ST-1.1.4)");
                    result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.4", "text": "Preferred", "ambiguous": false});
                    // unambiguous -> interpret...
                  } else {
                    // ST-1.1.4
                    // not acceptable
                    //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.1.4)");
                    result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.4", "text": "Not Acceptable", "ambiguous": true});
                    // ambiguous!
                  }
                }
              }


            } else if (numPlain == 2 && numWedge == 2 || numPlain == 2 && numHash == 2 || numWedge == 2 && numHash == 2) {
              var type1 = "", type2 = "";
              var sameType1 = [];
              var sameType2 = [];
              for (var j = 0; j < 4; j++) {
                if (!type1.length) {
                  type1 = bondTypes[j];
                } else if (!type2.length && type1 != bondTypes[j]) {
                  type2 = bondTypes[j]
                }
                if (type1 == bondTypes[j]) {
                  sameType1.push(j);
                } else if (type2 == bondTypes[j]) {
                  sameType2.push(j);
                }
              }

              if (!isAdjacent(sortedAngles, sameType1[0], sameType1[1])) {
                var angle1 = computeAngle(sortedAngles, sameType1[0], sameType2[0]);
                var angle2 = computeAngle(sortedAngles, sameType1[0], sameType2[1]);
                var angle3 = computeAngle(sortedAngles, sameType1[1], sameType2[0]);
                var angle4 = computeAngle(sortedAngles, sameType1[1], sameType2[1]);
                if (matchAngle(90.0, angle1) && matchAngle(90.0, angle2) && matchAngle(90.0, angle3) && matchAngle(90.0, angle4)) {
                  // ST-1.1.3
                  // acceptable
                  //JSChem.print("atom " + i + ": Acceptable (See ST-1.1.3)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.3", "text": "Acceptable", "ambiguous": false});
                  // unambiguous -> interpret...
                } else {
                  // ST-1.1.6
                  // not acceptable
                  //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.1.6)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.6", "text": "Not Acceptable", "ambiguous": false});
                  // unambiguous -> interpret...
                }

              } else {
                var largestGap = computeGap(sortedAngles);
                if (numPlain == 2) {
                  if (largestGap > 190.0) {
                    // ST-1.1.19
                    // wrong
                    //JSChem.print("atom " + i + ": Wrong (See ST-1.1.19)");
                    result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.19", "text": "Wrong", "ambiguous": true});
                    // ambiguous!
                  } else {
                    // ST-1.1.11
                    // not acceptable
                    //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.1.11)");
                    result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.11", "text": "Not Acceptable", "ambiguous": true});
                    // ambiguous!
                  }
                } else if (numPlain == 0) {
                  if (largestGap > 190.0) {
                    // ST-1.1.18
                    // wrong
                    //JSChem.print("atom " + i + ": Wrong (See ST-1.1.18)");
                    result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.18", "text": "Wrong", "ambiguous": true});
                    // ambiguous!
                  } else {
                    // ST-1.1.17
                    // wrong
                    //JSChem.print("atom " + i + ": Wrong (See ST-1.1.17)");
                    result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.17", "text": "Wrong", "ambiguous": true});
                    // ambiguous!
                  }
                }
              }

            } else if (numWedge == 3 || numHash == 3) {
              var type = numWedge == 3 ? "wedge" : "hash";
              var sameType = [];
              for (var j = 0; j < 4; j++) {
                if (bondTypes[j] == type) {
                  sameType.push(j);                  
                }
              }

              var angle1 = computeAngle(sortedAngles, sameType[0], sameType[1]);
              var angle2 = computeAngle(sortedAngles, sameType[1], sameType[2]);
              var angle3 = computeAngle(sortedAngles, sameType[0], sameType[2]);

              if (matchAngle(120.0, angle1) && matchAngle(120.0, angle2) && matchAngle(120.0, angle3)) {
                // ST-1.1.5
                // not acceptable
                //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.1.5)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.5", "text": "Not Acceptable", "ambiguous": false});
                // unambiguous -> interpret...
              } else {
                var bothOnSameSide = 0;
                for (var j = 0; j < 3; j++) {
                  var signs = [];
                  for (var k = 0; k < 3; k++) {
                    if (j != k) {
                      signs.push(getTriangleSign(bondVectors[sameType[j]], new vector2(0.0, 0.0), bondVectors[sameType[k]]));
                    }
                  }
                  if (signs[0] < 10.0 && signs[1] < 10.0 || signs[0] > -10.0 && signs[1] > -10.0) {
                    bothOnSameSide++;
                  }
                }
 
                if (bothOnSameSide > 1) {
                  // ST-1.1.9
                  // not acceptable
                  //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.1.9)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.9", "text": "Not Acceptable", "ambiguous": true});
                  // ambiguous!
                }


              }

            
            } else if (numPlain == 1 && (numWedge == 2 && numHash == 1 || numWedge == 1 && numHash == 2)) {
              // ST-1.1.7, ST-1.1.8 & ST-1.1.12
              // not acceptable
              //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.1.7, ST-1.1.8 & ST-1.1.12)");
              result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.[7,8,12]", "text": "Not Acceptable", "ambiguous": false});
              // unambiguous -> interpret...
            } else if (numWedge == 4 || numHash == 4) {
              // ST-1.1.20
              // wrong
              //JSChem.print("atom " + i + ": Wrong (See ST-1.1.20)");
              result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.1.20", "text": "Wrong", "ambiguous": true});
              // ambiguous!
            }

          //
          // ST-1.2 Tetrahedral configurations depicted with three explicit bonds
          // 
          } else if (numNbrs == 3) {
            if (numPlain == 2 && numWedge + numHash == 1) {
              var plain = [];
              var wedge_hash;
              for (var j = 0; j < 3; j++) {
                if (bondTypes[j] == "plain") {
                  plain.push(j);
                } else if (bondTypes[j] == "wedge") {
                  wedge_hash = j;
                } else if (bondTypes[j] == "hash") {
                  wedge_hash = j;              
                }
              }

              var plainAngle = computeAngle(sortedAngles, plain[0], plain[1]);
              
              if (plainAngle < 170.0) {
                var largestGap = computeGap(sortedAngles);
                if (largestGap > 180.0) {
                  var sign = getTriangleSign(new vector2(0.0, 0.0), bondVectors[wedge_hash], bondVectors[plain[0]]) *
                             getTriangleSign(new vector2(0.0, 0.0), bondVectors[wedge_hash], bondVectors[plain[1]]);
                  if (sign > 0.0) {
                    // ST-1.2.4
                    // no acceptable
                    //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.2.4)");
                    result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.4", "text": "Not Acceptable", "ambiguous": true});
                    // ambiguous! (unless in perspective ring, should be aboided though)
                  } else {
                    if (nbrIsInRing[plain[0]] && nbrIsInRing[plain[1]] && nbrIsInRing[wedge_hash]) {
                      // ST-1.2.3 exception (bridgehead atom)
                      // preferred
                      //JSChem.print("atom " + i + ": Preferred (See ST-1.2.3)");
                      result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.3", "text": "Preferred", "ambiguous": false});
                      // unambiguous -> interpret...
                    } else {
                      // ST-1.2.3
                      // not acceptable
                      //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.2.3)");
                      result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.3", "text": "Not Acceptable", "ambiguous": true});
                      // ambiguous!
                    }
                  }
                } else {
                  // ST-1.2.1
                  // preferred
                  //JSChem.print("atom " + i + ": Preferred (See ST-1.2.1)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.1", "text": "Preferred", "ambiguous": false});
                  // unambiguous -> interpret...
                }
              } else {
                // ST-1.2.12
                // wrong
                //JSChem.print("atom " + i + ": Wrong (See ST-1.2.12)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.12", "text": "Wrong", "ambiguous": true});
                // ambiguous!
              }
            } else if (numPlain == 1 && numWedge == 1 && numHash == 1) {
              var plain, wedge, hash;
              for (var j = 0; j < 3; j++) {
                if (bondTypes[j] == "plain") {
                  plain = j;
                } else if (bondTypes[j] == "wedge") {
                  wedge = j;
                } else if (bondTypes[j] == "hash") {
                  hash = j;              
                }
              }

              var wedgeHashAngle = computeAngle(sortedAngles, wedge, hash);

              if (matchAngle(30.0, wedgeHashAngle)) {
                // ST-1.2.2
                // acceptable
                //JSChem.print("atom " + i + ": Acceptable (See ST-1.2.2)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.2", "text": "Acceptable", "ambiguous": false});
                // unambiguous -> interpret...
              } else {
                // ST-1.2.10
                // wrong
                //JSChem.print("atom " + i + ": Wrong (See ST-1.2.10)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.10", "text": "Wrong", "ambiguous": true});
                // ambiguous!
              }

             
            } else if (numPlain == 1 && (numWedge == 2 || numHash == 2)) {

              // ST-1.2.5
              // not acceptable
              //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.2.5)");
              result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.5", "text": "Not Acceptable", "ambiguous": false});
              // unambiguous -> interpret...
             
            } else if (numWedge == 3 || numHash == 3) {
              var largestGap = computeGap(sortedAngles);
              if (largestGap < 180.0) {
                // ST-1.2.6
                // not acceptable
                //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.2.6)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.6", "text": "Not Acceptable", "ambiguous": false});
                // unambiguous -> interpret...
              } else { 
                // ST-1.2.9
                // not acceptable
                //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.2.9)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.9", "text": "Not Acceptable", "ambiguous": true});
                // ambiguous!
              }
 
            } else if (numWedge == 2 && numHash == 1 || numWedge == 1 && numHash == 2) {
              var wedge = [];
              var hash = [];
              for (var j = 0; j < 3; j++) {
                if (bondTypes[j] == "wedge") {
                  wedge.push(j);
                } else if (bondTypes[j] == "hash") {
                  hash.push(j);
                }
              }
 
              var largestGap = computeGap(sortedAngles);
              if (largestGap > 180.0) {
                var sameAngle;
                if (wedge.length == 2) {
                  sameAngle = computeAngle(sortedAngles, wedge[0], wedge[1]);
                } else if (hash.length == 2) {
                  sameAngle = computeAngle(sortedAngles, hash[0], hash[1]);
                }
                if (matchAngle(120.0, sameAngle)) {
                  // ST-1.2.7 
                  // not acceptable
                  //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.2.7)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.7", "text": "Not Acceptable", "ambiguous": false});
                  // unambiguous -> interpret...
                } else {
                  // ST-1.2.8
                  // not acceptable
                  //JSChem.print("atom " + i + ": Not Acceptable (See ST-1.2.8)");
                  result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.8", "text": "Not Acceptable", "ambiguous": true});
                  // ambiguous!
                }
              } else {
                // ST-1.2.11
                // wrong
                //JSChem.print("atom " + i + ": Wrong (See ST-1.2.11)");
                result.push({ "type": "tetrahedral", "index": i, "reference": "ST-1.2.11", "text": "Wrong", "ambiguous": true});
                // ambiguous!
              }

            
            }
          
          } // numNbrs == 3
        } // isTetrahedral
      } // foreach atom

      var numResults = result.length;
      for (var i = 0; i < numResults; i++) {
        var r = result[i];

        var text = r.text;
        if (text.indexOf("Not Acceptable") != -1) {
          text = "<b>" + text + "</b>";
        } else if (text.indexOf("Wrong") != -1) {
          text = "<b>" + text + "</b>";
        }

        if (r.type == "bond") {
          if (r.ambiguous) {
            JSChem.print("bond " + r.index + ": " + text + " (See " + r.reference + ") <font color=\"red\">Ambiguous!</font>");
          } else {
            JSChem.print("bond " + r.index + ": " + text + " (See " + r.reference + ")");
          }
        } else if (r.type == "tetrahedral") {
          if (r.ambiguous) {
            JSChem.print("atom " + r.index + ": " + text + " (See " + r.reference + ") <font color=\"red\">Ambiguous!</font>");
          } else {
            JSChem.print("atom " + r.index + ": " + text + " (See " + r.reference + ")");
          }
        }
      }

      return result;
    };

  };

  window.JSChem.plugin.register("stereo");

})();
