goog.provide('jchemhub.model.SingleBondUp');
goog.require('jchemhub.model.SingleBond');

/**
 * Class representing a Single Bond with Up sterochemistry
 * 
 * @param {jchemhub.model.Atom}
 *            source Atom at one end of bond.
 * @param {jchemhub.model.Atom}
 *            target Atom at other end of bond.
 * @constructor
 * @extends {jchemhub.model.SingleBond}
 */
jchemhub.model.SingleBondUp = function(source, target, opt_molecule) {
	jchemhub.model.SingleBond.call(this, source, target, opt_molecule);
}
goog.inherits(jchemhub.model.SingleBondUp, jchemhub.model.SingleBond);
/**
 * static value for order of this type of bond
 * 
 * @type{number}
 */
jchemhub.model.SingleBondUp.ORDER = 1;

/**
 * 
 * @return {jchemhub.model.SingleBondUp}
 */
jchemhub.model.SingleBondUp.prototype.clone = function() {
	return new jchemhub.model.SingleBondUp(this.source, this.target,
			this.molecule);
}