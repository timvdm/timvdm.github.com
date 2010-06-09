/**
 * <p>This page is meant for jchemhub developers to get started. It documents most
 * aspects of contributing code.</p>
 *
 * <h2 class="fixedFont">Getting Started</h2>
 * TODO: checking out source, running jslint, compiling with typechecking...
 *
 * <h2 class="fixedFont">Coding Style</h2>
 * <p>The jchemhub codng style uses a 4 space indentation. All modern editors can be 
 * configured to insert 4 spaces when pressing the tab button. All code is contained
 * in the jchemhub namespace and globals should be avoided at all cost. While the
 * javascript specification allows code to written without all ';' characters, it is
 * a good practise to put them in the correct place. Many of the style issues can be 
 * detected by running jslint. </p>
 * <p>There are also a few conventions regarding naming. Namespaces start with a lower
 * case letter. Classes start with a upper case letter. Methods start with lower case
 * letters. The other letters can be upper or lower case. Usually an upper case letter
 * is used when starting a new word. Some examples:</p>
 *
 * <pre class="code">
 * jchemhub.model.Atom
 * jchemhub.model.Atom.getNeighbors
 * jchemhub.symmetry.Symmetry.findSymmetryClasses
 * </pre>
 *
 * <p>Other parts of the code are often good examples. Another convention is the use
 * of interfaces. While javascript doesn't support interfaces, defining them makes 
 * it easier to document object requirements. The closure compiler also uses these
 * interfaces to do type checking. This makes the library more robust.</p>
 *
 * <h2 class="fixedFont">Documenting Code</h2>
 * <p>Documentation is an imortant part of a project. The documentation is generated 
 * from comments in the source code. This is done using a modified version of the 
 * <a href="http://code.google.com/p/jsdoc-toolkit/">jsdoc-toolkit</a>. Everything
 * works the same but some exclusive features are supported. For example, the @page
 * tag was added to enable custom pages to be added.</p>
 *
 * A typical code file could look like:
 * <pre class="code">
 * goog.provide('jchemhub.somenamespace.SomeClass')
 * goog.require('jchemhub.othernamespace.OtherClass')
 *
 * &#47;**
 *  * @namespace jchemhub
 *  * @name jchemhub
 *  *&#47;
 * &#47;**
 *  * @namespace jchemhub.somenamespace
 *  * @name jchemhub.somenamespace
 *  *&#47;
 *
 * &#47;**
 *  * 
 *  *
 *  *
 *  *       ... page with html markup ...
 *  *
 *  *
 *  * @page Some Page Name
 *  * @name needed, otherwise jsdoc will ignore the page...
 *  *&#47;
 * 
 * &#47;**
 *  * Description of the class goes here...
 *  *
 *  * Descriptions can be multi line...
 *  * &lt;b&gt;And html can be used&lt;/b&gt;
 *  *
 *  * @class 
 *  *&#47;
 * jchemhub.somenamespace.SomeClass = function() {
 *     &#47;**
 *      * Variable documentation
 *      *&#47;
 *     this.myVar = 0;
 * };
 *
 * &#47;**
 *  * Method description here...
 *  * 
 *  * Code example:
 *  * &lt;pre class="code"&gt;
 *  * var somevar = new jchemhub.somenamespace.SomeClass();
 *  * somevar.someMethod(3);
 *  * &lt;/pre&gt;
 *  *
 *  * @param value The value argument.
 *  * @return {Array.&lt;jchemhub.model.Atom&gt;} The atom for the specified value.
 *  *&#47;
 *  jchemhub.somenamespace.SomeClass.prototype.someMethod = function(&#47;**number*&#47;value) {
 *      ...
 *  };
 * </pre>
 *
 *
 *
 *
 * @page Developer Documentation
 * @name DeveloperDocumentation
 */
