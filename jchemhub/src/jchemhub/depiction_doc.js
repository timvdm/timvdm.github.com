/**
 *
 *
 * <h2 class="fixedFont">Introduction</h2>
 *
 * <p>JChemHub makes use of the Model View Controller design pattern. The models
 * are classes like {@link jchemhub.model.Molecule}, {@link jchemhub.model.Atom}
 * and {@link jchemhub.model.Bond}. The views are renderer classes to render
 * the various model classes. There are low level renderers for atoms, bonds 
 * and many other primitives. However, when using jchemhub as a library, only
 * the high level renders like {@link jchemhub.view.MoleculeRenderer} and 
 * {@link jchemhub.view.ReactionRenderer} are used directly.</p>
 *
 * <h2 class="fixedFont">Depicting a molecule</h2>
 *
 * <p>The google closure library is used for drawing the graphics. Depicting a
 * structure in a html page requires the following. A html div tag should be 
 * placed somewhere in the body where the depiction should go. The height and 
 * width can be specified using the style attribute. Make sure to give the div 
 * an id to retreive it.</p>
 * 
 * <pre class="code">
 * &lt;div id="reactionContainer" style="width: 640px; height: 480px"&gt;&lt;/div&gt;
 * </pre>
 *
 * <p>Create a function that will be called when the page has finished loading.
 * The relevant part of the function could look like:</p>
 * <pre class="code">
 * var graphicsContainer = goog.dom.getElement('graphicsContainer');
 * var graphics = goog.graphics.createGraphics(graphicsContainer.clientWidth, 
 *         graphicsContainer.clientHeight);
 * graphics.render(graphicsContainer);
 * var renderer = new jchemhub.view.MoleculeRenderer(null, graphics);
 * renderer.render(molecule);
 * </pre>
 * 
 * The example should be self explanatory. If not, the example below contains some comments:
 * 
 * <pre class="code">
 * &lt;!DOCTYPE html&gt;
 * &lt;html&gt;
 *     &lt;head&gt;
 *         &lt;title&gt;Depiction example&lt;/title&gt;
 *         &lt;script src="path/to/google/closure/goog/base.js"&gt;&lt;/script&gt;
 *         &lt;script src="path/to/jchemhub/deps.js"&gt;&lt;/script&gt;
 *         &lt;script&gt;
 *             goog.require('jchemhub.view.MoleculeRenderer');
 *             goog.require('goog.events.EventType'); // needed for load event
 *             goog.require('goog.dom'); // needed to get the graphicsContainer element
 *             goog.require('goog.graphics'); // needed for rendering
 *         &lt;/script&gt;
 *         &lt;script&gt;
 *             function initPage(){
 *                 var molecule = ...init molecule...;
 *                 // get the graphicsContainer
 *                 var graphicsContainer = goog.dom.getElement('graphicsContainer');
 *                 // create the goog.graphics object with required dimensions
 *                 var graphics = goog.graphics.createGraphics(graphicsContainer.clientWidth, 
 *                         graphicsContainer.clientHeight);
 *                 // set the graphics container
 *                 graphics.render(graphicsContainer);
 *                 // create the renderer
 *                 var renderer = new jchemhub.view.MoleculeRenderer(null, graphics);
 *
 *                 renderer.render(molecule);
 *             };
 *
 *             // call initPage when the page is loaded
 *             goog.events.listen(window, goog.events.EventType.LOAD, initPage);
 *         &lt;/script&gt;
 *     &lt;/head&gt;
 *     &lt;body&gt;
 *         &lt;div id="reactionContainer" style="width: 640px; height: 480px"&gt;&lt;/div&gt;
 *     &lt;/body&gt;
 * &lt;/html&gt;
 * </pre>
 *
 *
 * @page Depiction
 * @name Depiction
 */
