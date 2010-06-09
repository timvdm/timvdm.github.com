/**
 *
 *
 * <h2 class="fixedFont">Ring perception</h2>
 *
 * <p><b>Exhaustive Ring Perception</b></p>
 * <p>Exhaustive ring perception is done using the Hanser algorithm. The 
 * implementation is ported from the open source java cheminformatics 
 * library <a href="http://metamolecular.com/mx">MX from metamolecular</a>.
 * Reference:</p>
 *
 * <pre class="code">
 * A New Algorithm for Exhaustive Ring Perception in a Molecular Graph,
 * Th. Hanser, Ph. Jauffret, G. Kaufmann, J. Chem. Inf. Comput. Sci.,
 * 1996, 36(6), 1146-1152
 * <a href="http://pubs.acs.org/doi/abs/10.1021/ci960322f">link</a>
 * </pre>
 *
 * <p><b>Smallest Set of Smallest Rings (SSSR)</b></p>
 * <p>The SSSR algorithm uses path-included distance matrices to find the
 * Smallest Set of Smallest Rings. The first step of the algorithm computes
 * the distance matrix and two path matrices. The first path matrix contains
 * all shortest paths between atoms i and j in matrix element i,j. The second
 * path matrix contains the shortest+1 paths. From these matrices, the ring
 * set can be constructed. The algorithm is based on the Floyd-Warshall 
 * algorithm for computing the distance matrix. As a result, the runtime 
 * complexity is roughly O(n^3) where n is the number of atoms in the molecule.
 * Reference:</p>
 *
 * <pre class="code">
 * A robust method for searching the smallest set of smallest rings with a 
 * path-included distance matrix, Chang Lee, PNAS, Oct. 13, 2009, vol. 16,
 * no. 41, pages 17355-17358
 * <a href="http://www.pnas.org/content/106/41/17355.full">paper (open access)</a>
 * <a href="http://www.pnas.org/content/106/41/17355/suppl/DCSupplemental">supplementary information (pseudo code for the algorithm</a>
 * </pre>
 *
 * <h2 class="fixedFont">Symmetry Classes</h2>
 * <p>Symmetry classes are used to represent graph symmetry in a molecule's topology.
 * In graph theory, these symmetry classes are often called colors. Every node 
 * (i.e. atom) gets a symmetry class or color assigned to it and symmetric atoms 
 * will get the same symmetry class. Such atoms are also called isospectral points.
 * Below are two molecules annotated with symmetry classes. The left molecule has 
 * one symmetry axis, the right has two.</p>
 * <center><img src="symmetryclasses.png"></center>
 * <p>The algorithm for assigning symmetry classes is bases on Morgan's Extended 
 * Connectivity algorithm. The first step is to assign graph invariants to each atom.
 * A graph invariant can be anything that doesn't rely on atom ordering and 
 * discriminates between different atoms (i.e. atoms with the same chemical and 
 * topological properties). An example of a graph invariant is an atom's degree.
 * However, using atom degree as initial values, false isospectral atoms will be 
 * obtained. The discriminating power of the graph invariant can be improved by 
 * using the graph theoretical distance. This is the sum of the distances from an 
 * atom to all other atoms in the molecule. Using this as graph invariant to encode
 * topology, acceptable results are obtained for most chemical graphs. Other 
 * topological and chemical information such as atomic number, ring membership
 * are also included.</p>
 * <p>Once the graph invariants are assigned, an iterative process begins to
 * find the final values. In each iteration, the atom values are replaced by
 * the sum of the atom's neighbor values. The iteration goes on until the 
 * number of unique symmetry classes remains the same after an iteration. Next,
 * The symmetry classes are renumbered from 1 to x where x is the number of
 * unique symmetry classes in the molecule.
 *
 * @page Algorithms
 * @name Algorithms
 */
