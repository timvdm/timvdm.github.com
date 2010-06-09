goog.provide('jchemhub.model.TripleBond');
goog.require('jchemhub.model.Bond');

/**
 * Class representing a Triple Bond
 * 
 * @param {jchemhub.model.Atom}
 *            source Atom at one end of bond.
 * @param {jchemhub.model.Atom}
 *            target Atom at other end of bond.
 * @constructor
 * @extends {jchemhub.model.Bond}
 */
jchemhub.model.TripleBond = function(source, target, opt_molecule) {
	jchemhub.model.Bond.call(this, source, target, opt_molecule);
}
goog.inherits(jchemhub.model.TripleBond, jchemhub.model.Bond);

/**
 * static value for order of this type of bond
 * 
 * @type{number}
 */
jchemhub.model.TripleBond.ORDER = 3;

jchemhub.model.TripleBond.prototype.clone = function() {
	return new jchemhub.model.TripleBond(this.source, this.target,
			this.molecule);
}