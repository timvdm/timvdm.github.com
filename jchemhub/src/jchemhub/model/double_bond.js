goog.provide('jchemhub.model.DoubleBond');
goog.require('jchemhub.model.Bond');

/**
 * Class representing a Double Bond
 * 
 * @param {jchemhub.model.Atom}
 *            source Atom at one end of bond.
 * @param {jchemhub.model.Atom}
 *            target Atom at other end of bond.
 * @constructor
 * @extends {jchemhub.model.Bond}
 */
jchemhub.model.DoubleBond = function(source, target, opt_molecule) {
	jchemhub.model.Bond.call(this, source, target, opt_molecule);
}
goog.inherits(jchemhub.model.DoubleBond, jchemhub.model.Bond);

/**
 * static value for order of this type of bond
 * 
 * @type{number}
 */
jchemhub.model.DoubleBond.ORDER = 2;

/**
 * 
 * @param {jchemhub.model.DoubleBond} bond Bond to clone.
 * @return {jchemhub.model.DoubleBond}
 */
jchemhub.model.DoubleBond.prototype.clone = function() {
	return new jchemhub.model.DoubleBond(this.source, this.target,
			this.molecule);
}
