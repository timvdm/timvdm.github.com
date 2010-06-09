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

goog.provide('jchemhub.query.DFSMapper');

goog.require('jchemhub.ring.SSSRFinder');
goog.require('goog.structs.Map');

(function() {

    function simpleSort(a, b) {
        return a - b
    }

    /**
     * This class implements the {@link jchemhub.query.IMapper} interface and can
     * be used for substructure searching. The algorithm compares each query atom
     * with the molecule atoms in a depth-first order. As long as the atom pairs
     * match, the path is extended. When the path contains the same number of 
     * elements as there are atoms in the query, a match is found. See the
     * <a href="../Substructure Search.html">Substructure Search</a> page for more
     * details.
     * @class Class implementing the jchemhub.query.IMapper interface.
     * @implements {jchemhub.query.IMapper}
     */
    jchemhub.query.DFSMapper = function(query) {

        var Type = { MapFirst: 0, MapUnique: 1, MapAll: 2 };

        this.query = query;

        /**
         * Select a start atom. todo: select the most unique query atom,
         * this avoids checking many paths that will never lead to a match. 
         */
        function getStartAtom(query) {
            return query.getAtom(0);
        }

        /**
         * State to represent the current mapped state.
         */
        function State(type, query, queried) {
            this.type = type; // MapFirst, MapUnique, MapAll
            this.query = query; // the query
            this.queried = queried; // the queried molecule
            this.sssr = jchemhub.ring.findSSSR(queried); // sssr rings
            this.queryPath = []; // the path in the query
            this.queriedPath = []; // the path in the queried molecule
            this.candidates = []; // the candidates to check
        }

        /**
         * Check if the current state is a full mapping of the query.
         */
        function checkForMap(state, maps) {
            // store the mapping if all atoms are mapped
            if (state.queryPath.length === state.query.countAtoms()) {
                // create the map
                var map = new goog.structs.Map();
                for (var k = 0, kl = state.queryPath.length; k < kl; k++) {
                    map.set(state.query.indexOfAtom(state.queryPath[k]), state.queried.indexOfAtom(state.queriedPath[k]));
                }
                if (state.type === Type.MapUnique) {
                    var values = map.getValues();
                    values.sort(simpleSort);
                    var isUnique = true;
                    for (k = 0, kl = maps.length; k < kl; k++) {
                        var kValues = maps[k].getValues();
                        kValues.sort(simpleSort);
                        if (goog.array.equals(values, kValues)) {
                            isUnique = false;
                        }
                    }
                    if (isUnique) {
                        maps.push(map);
                    }
                } else {
                    maps.push(map);
                }

                if (state.type === Type.MapFirst) {
                    return;
                }
            }
        }

        /**
         * Match the candidate atoms and bonds.
         */
        function matchCandidate(state, queryAtom, queriedAtom, queryNbr, queriedNbr, maps) {
            // make sure the neighbor atom isn't in the paths already
            if (goog.array.indexOf(state.queryPath, queryNbr) !== -1) {
                return false;
            }
            if (goog.array.indexOf(state.queriedPath, queriedNbr) !== -1) {
                return false;
            }

            // check if the atoms match
            if (!queryNbr.matches(queriedNbr, state.queried, state.sssr)) {
                return false;
            }

            var queryBond = state.query.findBond(queryAtom, queryNbr);
            var queriedBond = state.queried.findBond(queriedAtom, queriedNbr);

            // check if the bonds match
            if (!queryBond.matches(queriedBond, state.queried, state.sssr)) {
                return false;
            }

            // add the neighbors to the paths
            state.queryPath.push(queryNbr);
            state.queriedPath.push(queriedNbr);

            // check if this is a full match
            checkForMap(state, maps);

            return true;
        }

        /**
         * The depth-first isomorphism algorithm.
         */
        function DFS(state, queryAtom, queriedAtom, maps) {
            var queryNbrs = queryAtom.getNeighbors();
            var queriedNbrs = queriedAtom.getNeighbors();

            // load the possible candidates
            var candidates = [];
            for (var i = 0, li = queryNbrs.length; i < li; i++) {
                var queryNbr = queryNbrs[i];
                for (var j = 0, lj = queriedNbrs.length; j < lj; j++) {
                    var queriedNbr = queriedNbrs[j];
                    candidates.push({ queryAtom: queryAtom, queryNbr: queryNbr, queriedAtom: queriedAtom, queriedNbr: queriedNbr });
                }
            }
            // save the candidates for later (used to explore branches)
            if (state.candidates.length) {
                state.candidates.push(state.candidates[state.candidates.length-1].concat(candidates));
            } else {
                state.candidates.push(candidates);
            }

            // do the mapping by checking the candidates
            while (state.candidates[state.candidates.length-1].length) {
                var candidate = state.candidates[state.candidates.length-1].pop();
                if (matchCandidate(state, candidate.queryAtom, candidate.queriedAtom, candidate.queryNbr, candidate.queriedNbr, maps)) {
                    DFS(state, candidate.queryNbr, candidate.queriedNbr, maps);

                    // backtrack
                    state.queryPath.pop();
                    state.queriedPath.pop();
                    state.candidates.pop();
                }
            }
        }

        /**
         * Helper function for mapAllCallback, mapUniqueCallback and MapFirstCallback
         */
        function mapNext(i, type, query, queryAtom, queried, maps, callback) {
            var state = new State(type, query, queried);
            var queriedAtom = queried.getAtom(i);

            if (queryAtom.matches(queriedAtom)) {
                state.queryPath.push(queryAtom);
                state.queriedPath.push(queriedAtom);
                DFS(state, queryAtom, queriedAtom, maps);
            }

            i++;
            if (i < queried.countAtoms()) {
                var nextBitOfWork = function() { mapNext(i, type, query, queryAtom, queried, maps, callback); };
                setTimeout(nextBitOfWork, 0);
            } else {
                callback(maps);
            }
        }

        /**
         * Get all mappings of the query on the queried molecule.
         * @param {jchemhub.model.Molecule} queried The queried molecule.
         * @return {Array.<goog.structs.Set>} The mappings
         */
        this.mapAll = function(queried) {
            var maps = [];
            var queryAtom = getStartAtom(this.query);
            for (var i = 0, li = queried.countAtoms(); i < li; i++) {
                var state = new State(Type.MapAll, this.query, queried);
                var queriedAtom = queried.getAtom(i);
                if (!queryAtom.matches(queriedAtom)) {
                    continue;
                }

                if (this.query.countAtoms() > 1) {
                    state.queryPath.push(queryAtom);
                    state.queriedPath.push(queriedAtom);
                    DFS(state, queryAtom, queriedAtom, maps);
                } else {
                    var map = new goog.structs.Map();
                    map.set(state.query.indexOfAtom(queryAtom), state.queried.indexOfAtom(queriedAtom));
                    maps.push(map);
                }
            }

            return maps;
        };

        /**
         * Get all mappings of the query on the queried molecule. The specified
         * callback function will be called with the found maps as argument.
         * Unlike mapAll, this function regulary gives control back to the
         * browser, preventing the GUI to lock up.
         * @param {jchemhub.model.Molecule} queried The queried molecule.
         * @param {Function} callback The callback function to report results.
         */
        this.mapAllCallback = function(queried, callback) {
            var maps = [];
            var queryAtom = getStartAtom(this.query);
            var i = 0;
            mapNext(i, Type.MapAll, this.query, queryAtom, queried, maps, callback);
        };

        /**
         * Get all unique mappings of the query on the queried molecule.
         * @param {jchemhub.model.Molecule} queried The queried molecule.
         * @return {Array.<goog.structs.Set>} The unique mappings
         */
        this.mapUnique = function(queried) {
            var maps = [];
            var queryAtom = getStartAtom(this.query);
            for (var i = 0, li = queried.countAtoms(); i < li; i++) {
                var state = new State(Type.MapUnique, this.query, queried);
                var queriedAtom = queried.getAtom(i);
                if (!queryAtom.matches(queriedAtom)) {
                    continue;
                }
                if (this.query.countAtoms() > 1) {
                    state.queryPath.push(queryAtom);
                    state.queriedPath.push(queriedAtom);
                    DFS(state, queryAtom, queriedAtom, maps);
                } else {
                    var map = new goog.structs.Map();
                    map.set(state.query.indexOfAtom(queryAtom), state.queried.indexOfAtom(queriedAtom));
                    maps.push(map);
                }
            }

            return maps;
        };

        /**
         * Get all unique mappings of the query on the queried molecule. The
         * specified callback function will be called with the found maps as
         * argument. Unlike mapUnique, this function regulary gives control
         * back to the browser preventing the GUI to lock up.
         * @param {jchemhub.model.Molecule} queried The queried molecule.
         * @param {Function} callback The callback function to report results.
         */
        this.mapUniqueCallback = function(queried, callback) {
            var maps = [];
            var queryAtom = getStartAtom(this.query);
            var i = 0;
            mapNext(i, Type.MapUnique, this.query, queryAtom, queried, maps, callback);
        };


        /**
         * Get the first mappings of the query on the queried molecule.
         * @param {jchemhub.model.Molecule} queried The queried molecule.
         * @return {goog.structs.Set} The mapping
         */
        this.mapFirst = function(queried) {
            var maps = [];
            var queryAtom = getStartAtom(this.query);
            for (var i = 0, li = queried.countAtoms(); i < li; i++) {
                var state = new State(Type.MapFirst, this.query, queried);
                var queriedAtom = queried.getAtom(i);
                if (!queryAtom.matches(queriedAtom)) {
                    continue;
                }

                if (this.query.countAtoms() > 1) {
                    state.queryPath.push(queryAtom);
                    state.queriedPath.push(queriedAtom);
                    DFS(state, queryAtom, queriedAtom, maps);
                } else {
                    var map = new goog.structs.Map();
                    map.set(state.query.indexOfAtom(queryAtom), state.queried.indexOfAtom(queriedAtom));
                    return map;
                }

                if (maps.length) {
                    return maps[0];
                }
            }

            return new goog.structs.Map();
        };

    };


}());
